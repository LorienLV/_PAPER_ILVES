#pragma once

#include "graph.h"

#include "gromacs/topology/idef.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/real.h"
#include "gromacs/domdec/domdec.h"

/**
 * The Molecule structure contains all the information about the bonds and atoms
 * of the molecule required by ILVES constraint solver. A Molecule structure can
 * contain several molecules, with out any connection between them.
 *
 */

namespace gmx {

class Molecule {
public:
    struct Atoms {
        int num;       // The number of atoms in the molecule. Only atoms with
                       // constraints.
        Graph graph;   // the atomic graph of the molecule. Only atoms with
                       // constraints.

        const real *invmass;   // invmass[i] is 1/mass of the ith atom
    } atoms;

    struct Bonds {
        int num;       // The number of bonds in the molecule.
        Graph graph;   // The bond graph of the molecule.

        // GROMACS atom indices of the two atoms of each bond. Example:
        // Bond 0 connects atom1[0] and atom2[0].
        std::vector<int, AlignedAllocator<int>> atom1;
        std::vector<int, AlignedAllocator<int>> atom2;

        // Local atom indices of the two atoms of each bond.
        // The local atom indices correspond to the numbering of
        // atoms_data.graph. Example: Bond 0 connects latom1[0] and latom2[0].
        std::vector<int, AlignedAllocator<int>> latom1;
        std::vector<int, AlignedAllocator<int>> latom2;

        // std::vector<real, AlignedAllocator<real>> sigmaA;
        // std::vector<real, AlignedAllocator<real>> sigmaB;

        // The bond length squared of each bond.
        std::vector<real, AlignedAllocator<real>> sigma2;
        std::vector<real, AlignedAllocator<real>> invsigma2;
    } bonds;

    /**
     * Default constructor is deleted to avoid errors.
     */
    Molecule() = delete;

    /**
     * Constructs a Molecule that can be used by the ILVES constraint solver.
     *
     * @param idef The interaction definitions provided by GROMACS.
     * @param cr MPI communication structure.
     * @param invmass A pointer to the array of inverse mass of each atom.
     * The Molecule keeps a pointer to this array, so it must be kept alive.
     */
    Molecule(const InteractionDefinitions &idef, const t_commrec *cr, const real *invmass);

    /**
     * Renumber the data of the Bonds structure given a permutation.
     * The permutation is given as in MATLAB. Example:
     *  p = [2, 1, 0] Means that
     *  Old position 2 is now position 0
     *  Old position 1 is now position 1
     *  Old position 0 is now position 2
     *
     * @param perm The permutation in MATLAB format.
     * @param renumber_graph If true, the bond graph is also renumbered. The
     * bond graph is not renumbered otherwise.
     */
    void renumber_bonds(const std::vector<int> &perm, bool renumber_graph);
};
};   // namespace gmx
