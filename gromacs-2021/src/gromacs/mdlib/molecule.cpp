#include "molecule.h"

#include <vector>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/real.h"
#include "gromacs/topology/ifunc.h"
#include <algorithm>

#include <limits>

namespace gmx {

Molecule::Molecule(const InteractionDefinitions &idef,
                   const t_commrec *const cr,
                   const real *const _invmass) {

    atoms.invmass = _invmass;

    // The number of atoms with constraints.
    atoms.num = 0;

    // The number of local (rank) constraints.
    bonds.num = idef.il[F_CONSTR].size() / 3;

    bonds.atom1.resize(bonds.num);
    bonds.atom2.resize(bonds.num);

    bonds.latom1.resize(bonds.num);
    bonds.latom2.resize(bonds.num);

    bonds.sigma2.resize(bonds.num);
    bonds.invsigma2.resize(bonds.num);

    /*
     * OpenMP global variables
     */
    // The bond that connects each pair of atoms.
    std::vector<int> bonds_of_atom;

    int max_latom_id = std::numeric_limits<int>::min();
    int min_latom_id = std::numeric_limits<int>::max();

    #pragma omp parallel
    {
        const auto &iatoms = idef.il[F_CONSTR].iatoms;

        // Separate the bond lengths into two arrays to be SIMD friendly
        // and precompute sigma2.
        #pragma omp for reduction(max:max_latom_id) reduction(min:min_latom_id)
        for (int bond = 0; bond < bonds.num; ++bond) {
            const int type = iatoms[bond * 3 + 0];
            // Global atom indices.
            const int a = iatoms[bond * 3 + 1];
            const int b = iatoms[bond * 3 + 2];

            bonds.atom1[bond] = a;
            bonds.atom2[bond] = b;

            // mol->sigmaA[bond] = idef.iparams[type].constr.dA;
            // mol->sigmaB[bond] = idef.iparams[type].constr.dB;
            bonds.sigma2[bond] = square(idef.iparams[type].constr.dA);
            bonds.invsigma2[bond] = 1.0 / bonds.sigma2[bond];

            int la = a;
            int lb = b;

            // If domain decomposition, get the topology indexes of the atoms.
            if (DOMAINDECOMP(cr)) {
                la = cr->dd->globalAtomIndices[a];
                lb = cr->dd->globalAtomIndices[b];
            }

            bonds.latom1[bond] = la;
            bonds.latom2[bond] = lb;

            const int min_la_lb = std::min(la, lb);
            const int max_la_lb = std::max(la, lb);

            if (min_la_lb < min_latom_id) {
                min_latom_id = min_la_lb;
            }
            if (max_la_lb > max_latom_id) {
                max_latom_id = max_la_lb;
            }
        }
        // Implicit wait.

        #pragma omp master
        {
            // Resize atom.graph.xadj.
            atoms.graph.xadj.resize(max_latom_id - min_latom_id + 2, -2);

            // Now xadj[i] > -2 if i is in latom1 or latom2.
            for (int bond = 0; bond < bonds.num; ++bond) {
                const int la = bonds.latom1[bond];
                const int lb = bonds.latom2[bond];

                atoms.graph.xadj[la - min_latom_id] = -1;
                atoms.graph.xadj[lb - min_latom_id] = -1;
            }

            // Set the local id, such that the ids of the local atoms are
            // contiguous.
            atoms.num = 0;
            for (int idx = 0; idx < max_latom_id - min_latom_id + 1; ++idx) {
                if (atoms.graph.xadj[idx] == -1) {
                    atoms.graph.xadj[idx] = atoms.num;
                    ++atoms.num;
                }
            }
        }

        #pragma omp barrier

        // Apply the new numbering to the local atom ids.
        #pragma omp for
        for (int bond = 0; bond < bonds.num; ++bond) {
            bonds.latom1[bond] = atoms.graph.xadj[bonds.latom1[bond] - min_latom_id];
            bonds.latom2[bond] = atoms.graph.xadj[bonds.latom2[bond] - min_latom_id];
        }
        // Implicit wait.

        // Reinitialize atom.graph.xadj. See later why it is initialized with
        // 1s.
        #pragma omp for
        for (int atom = 0; atom < atoms.num; ++atom) {
            atoms.graph.xadj[atom + 1] = 1;
        }
        // Implicit wait.

        // Construct the atomic graph.

        // Compute xadj and adj.
        #pragma omp master
        {
            atoms.graph.nnodes = atoms.num;
            atoms.graph.xadj.resize(atoms.num + 1);

            // In the following loop, xadj[i + 1] will be the number of
            // connections of the ith atom. xadj[0] is not used. Array
            // initialized with 1s since each atom is connected to itself.
            for (int bond = 0; bond < bonds.num; ++bond) {
                const int la = bonds.latom1[bond];
                const int lb = bonds.latom2[bond];

                // Atom la is connected to lb (+1 number of connections).
                ++atoms.graph.xadj[la + 1];
                // Same for lb.
                ++atoms.graph.xadj[lb + 1];
            }

            // Process the current xadj to be the final xadj.
            atoms.graph.xadj[0] = 0;
            for (int atom = 0; atom < atoms.num; ++atom) {
                atoms.graph.xadj[atom + 1] = atoms.graph.xadj[atom] +
                                             atoms.graph.xadj[atom + 1];
            }

            atoms.graph.adj.resize(atoms.graph.xadj[atoms.num]);

            // The index of the next element to be inserted in the adj of each
            // atom.
            std::vector<int> adj_idx(atoms.num, 0);

            bonds_of_atom.resize(atoms.graph.xadj[atoms.num]);

            // Fill adj.
            for (int bond = 0; bond < bonds.num; ++bond) {
                const int la = bonds.latom1[bond];
                const int lb = bonds.latom2[bond];

                // Atom la is connected to lb.
                atoms.graph.adj[atoms.graph.xadj[la] + adj_idx[la]] = lb;
                // Same for lb.
                atoms.graph.adj[atoms.graph.xadj[lb] + adj_idx[lb]] = la;

                // The bond that connects la and lb is bond.
                bonds_of_atom[atoms.graph.xadj[la] + adj_idx[la]] = bond;
                bonds_of_atom[atoms.graph.xadj[lb] + adj_idx[lb]] = bond;

                ++adj_idx[la];
                ++adj_idx[lb];
            }
        }

        #pragma omp barrier // Wait for adj to be computed.

        // Sort adj.
        #pragma omp for nowait
        for (int atom = 0; atom < atoms.num; ++atom) {
            // Atoms are connected to themselves.
            atoms.graph.adj[atoms.graph.xadj[atom + 1] - 1] = atom;
            // But this connection is not a bond.
            // bonds_of_atom[atoms.graph.xadj[atom + 1] - 1] = -1;

            // The atom ids in adj are sorted.
            std::sort(atoms.graph.adj.begin() + atoms.graph.xadj[atom],
                      atoms.graph.adj.begin() + atoms.graph.xadj[atom + 1]);
        }

        // Construct the bond graph.

        // Resize xadj and adj. Compute xadj.
        #pragma omp master
        {
            bonds.graph.nnodes = bonds.num;

            bonds.graph.xadj.resize(bonds.num + 1);
            bonds.graph.xadj[0] = 0;

            for (int bond = 0; bond < bonds.num; ++bond) {
                const int la = bonds.latom1[bond];
                const int lb = bonds.latom2[bond];

                // The bond is connected to all bonds connected to atom la and
                // lb plus itself. -2 since both la and lb are connected between
                // them (this is the bond; counted in nconn with +1) and to
                // themselves (this is not a bond).
                const int nbonds_la = atoms.graph.xadj[la + 1] -
                                      atoms.graph.xadj[la] - 2;
                const int nbonds_lb = atoms.graph.xadj[lb + 1] -
                                      atoms.graph.xadj[lb] - 2;
                const int nconn = nbonds_la + nbonds_lb + 1;

                bonds.graph.xadj[bond + 1] = bonds.graph.xadj[bond] + nconn;
            }

            bonds.graph.adj.resize(bonds.graph.xadj.back());
        }

        #pragma omp barrier // Wait for xadj and adj to be resized.

        // Fill adj.
        #pragma omp for
        for (int bond = 0; bond < bonds.num; ++bond) {
            const int la = bonds.latom1[bond];
            const int lb = bonds.latom2[bond];

            // The bond is connected to all bonds connected to atom la and lb
            // plus itself.
            int adj_idx = bonds.graph.xadj[bond];
            // The bond is connected to itself.
            bonds.graph.adj[adj_idx] = bond;
            ++adj_idx;

            auto bonds_of_atom_to_adj = [&](const int atom) {
                for (int k = atoms.graph.xadj[atom];
                     k < atoms.graph.xadj[atom + 1] - 1;
                     ++k) {
                    const int neigh = bonds_of_atom[k];

                    // The bond is connected to itself only once.
                    if (neigh == bond) {
                        continue;
                    }

                    bonds.graph.adj[adj_idx] = neigh;

                    ++adj_idx;
                }
            };

            // The bond is connected to all bonds connected to atom la and lb.
            bonds_of_atom_to_adj(la);
            bonds_of_atom_to_adj(lb);

            // The bonds ids are sorted.
            std::sort(bonds.graph.adj.begin() + bonds.graph.xadj[bond],
                      bonds.graph.adj.begin() + bonds.graph.xadj[bond + 1]);
        }
    }
}

void Molecule::renumber_bonds(const std::vector<int> &perm,
                              const bool renumber_graph) {
    Molecule::Bonds new_bonds;

    // Copy old vectors.
    new_bonds.atom1 = bonds.atom1;
    new_bonds.atom2 = bonds.atom2;

    new_bonds.latom1 = bonds.latom1;
    new_bonds.latom2 = bonds.latom2;

    new_bonds.sigma2 = bonds.sigma2;
    new_bonds.invsigma2 = bonds.invsigma2;

    // Apply permutation to the vectors.
    for (int bond = 0; bond < bonds.num; bond++) {
        bonds.atom1[bond] = new_bonds.atom1[perm[bond]];
        bonds.atom2[bond] = new_bonds.atom2[perm[bond]];

        bonds.latom1[bond] = new_bonds.latom1[perm[bond]];
        bonds.latom2[bond] = new_bonds.latom2[perm[bond]];

        bonds.sigma2[bond] = new_bonds.sigma2[perm[bond]];
        bonds.invsigma2[bond] = new_bonds.invsigma2[perm[bond]];
    }

    // Renumber the vertices of the bond graph.
    if (renumber_graph) {
        std::vector<int> iperm;
        bonds.graph.renumber_vertices(perm, iperm);
    }
}

};   // namespace gmx