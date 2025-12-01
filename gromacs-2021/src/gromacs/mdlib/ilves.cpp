#include "ilves.h"

#include <functional>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_simd.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx {

Ilves::Ilves(const t_commrec *const cr,
             const InteractionDefinitions &idef,
             const real *invmass,
             const int threads,
             const bool upper_tri,
             std::vector<std::map<std::string, std::chrono::duration<double, std::milli>>>
                 &durations)
    : cr(cr), nthreads(threads), current_lagr(threads + 1), durations(durations) {

#if PRINT_FUNCTION_DURATIONS
    durations.resize(nthreads);
#endif

#if PRINT_FUNCTION_DURATIONS
    auto start = std::chrono::high_resolution_clock::now();
#endif
    mol = std::make_unique<Molecule>(idef, cr, invmass);
#if PRINT_FUNCTION_DURATIONS
    auto end = std::chrono::high_resolution_clock::now();
    durations[0]["Init Molecule"] = end - start;
#endif

#if PRINT_FUNCTION_DURATIONS
    start = std::chrono::high_resolution_clock::now();
#endif

    for (int d = 0; d < DIM; ++d) {
        x_ab[d].resize(mol->bonds.num);
    }

    const int nparts = nthreads + 1;

    std::vector<int> partitioning(nparts + 1);

    std::vector<int> bonds_perm_kpart;

    /*
     * Partitioning of the bonds between threads.
     */
    if (nthreads > 1) {
        // Partition id of each bond.
        std::vector<int> bonds_part_id;
        // Number of bonds per partition.
        std::vector<int> nbonds_per_part(nparts, 0);

        // Heuristic to determine the maximum allowed size of a submolecules or
        // blocks to keep load balance.
        const int block_max_size = std::min(10, mol->bonds.num / nthreads);
        const bool is_disjoint_mol = disjoint_mol(block_max_size);

        if (is_disjoint_mol) {
            // We can do a k-way partitioning of the **bond graph** without
            // any cut.
            bonds_part_id = mol->bonds.graph.kway_partition_disjoint(nthreads);
            for (int bond = 0; bond < mol->bonds.num; ++bond) {
                ++nbonds_per_part[bonds_part_id[bond]];
            }
        }
        else {
            // Partitioning of the **atomic** graph.
            const auto atoms_part_id = mol->atoms.graph.kway_partition(nthreads);

            bonds_part_id.resize(mol->bonds.num);

            // Determine which edges (bonds) has been cut.
            for (int bond = 0; bond < mol->bonds.graph.num_nodes(); ++bond) {
                // Get the pair of local atoms of the bond.
                const int la = mol->bonds.latom1[bond];
                const int lb = mol->bonds.latom2[bond];

                // If both atoms are not in the same partition, this a
                // separation bond.
                const int part = (atoms_part_id[la] != atoms_part_id[lb])
                                     ? nthreads
                                     : atoms_part_id[la];

                bonds_part_id[bond] = part;
                ++nbonds_per_part[part];
            }
        }

        // We can now construct the partitioning of the arrowhead matrix
        partitioning[0] = 0;
        for (int part = 0; part < nparts; ++part) {
            partitioning[part + 1] = partitioning[part] + nbonds_per_part[part];
        }

        // Permute the numbering of the bonds to have consecutive bonds in each
        // partition.

        // Compute the first bond idx of each partition.
        std::vector<int> bonds_part_idx(nparts, 0);
        for (int p = 1; p < nparts; ++p) {
            bonds_part_idx[p] = bonds_part_idx[p - 1] + nbonds_per_part[p - 1];
        }

        bonds_perm_kpart.resize(mol->bonds.num);
        std::vector<int> bonds_iperm_kpart(mol->bonds.num);

        // First compute the inverse permutation.
        for (int bond = 0; bond < mol->bonds.num; ++bond) {
            const int part = bonds_part_id[bond];

            bonds_iperm_kpart[bond] = bonds_part_idx[part];
            ++bonds_part_idx[part];
        }

        // Now we can compute the final permutation.
        for (int bond = 0; bond < mol->bonds.num; ++bond) {
            bonds_perm_kpart[bonds_iperm_kpart[bond]] = bond;
        }

        mol->bonds.graph.renumber_vertices(bonds_perm_kpart, bonds_iperm_kpart);
    }
    else {
        partitioning[0] = 0;
        partitioning[1] = mol->bonds.num;
        partitioning[2] = partitioning[1];
    }

#if PRINT_FUNCTION_DURATIONS
    end = std::chrono::high_resolution_clock::now();
    durations[0]["Init Partitioning"] = end - start;
#endif

#if PRINT_FUNCTION_DURATIONS
    start = std::chrono::high_resolution_clock::now();
#endif

    std::vector<int> schur_solver_perm;
    schur_solver = std::make_unique<SchurLinearSolver>(mol->bonds.graph,
                                                       upper_tri,
                                                       partitioning,
                                                       schur_solver_perm);
#if PRINT_FUNCTION_DURATIONS
    end = std::chrono::high_resolution_clock::now();
    durations[0]["Init Schur solver"] = end - start;
#endif

#if PRINT_FUNCTION_DURATIONS
    start = std::chrono::high_resolution_clock::now();
#endif

    // Construct the final permutation.
    std::vector<int> bonds_perm;
    if (!bonds_perm_kpart.empty()) {
        bonds_perm.resize(mol->bonds.num);
        for (int bond = 0; bond < mol->bonds.num; ++bond) {
            bonds_perm[bond] = bonds_perm_kpart[schur_solver_perm[bond]];
        }
    }
    else {
        bonds_perm = std::move(schur_solver_perm);
    }
    // WARNING: We do not need to use mol->bonds.graph again, so save time by
    // avoiding renumbering it. Also, the bond graph has already been renumbered
    // after computing the partition, so applying the permutation again would
    // be wrong.
    mol->renumber_bonds(bonds_perm, false);

    make_weights();

    // Resize the current_lagr vectors.
    // It needs to be as big as rhs in order to swap them at some point.
    for (int p = 0; p < schur_solver->part_data.size(); ++p) {
        current_lagr[p].resize(schur_solver->part_data[p].rhs.size());

        // Print the size of the partition
        // std::cout << "Partition " << p << " has "
        //           << schur_solver->part_data[p].rhs.size() << " rows.\n";
    }

#if PRINT_FUNCTION_DURATIONS
    end = std::chrono::high_resolution_clock::now();
    durations[0]["Init Permutations"] = end - start;
#endif

#if 0
    if (DOMAINDECOMP(cr)) {
        if (cr->dd->constraints) {
            dd_get_constraint_range(cr->dd,
                                    &non_local_atoms_range.first,
                                    &non_local_atoms_range.second);
        }
        else {
            int home_atoms = dd_numHomeAtoms(*cr->dd);
            non_local_atoms_range = std::make_pair(home_atoms, home_atoms);
        }
    }
    else {
        // No non-local atoms.
        non_local_atoms_range = std::make_pair(mol->gm, mol->gm);
    }
#endif
}

namespace {

/**
 * Subtraction of two rvec taking the PBC into account if pbc is not null.
 *
 * @param pbc The PBC (container) information. Null if there is no PBC
 * information.
 * @param vin1 First source rvec.
 * @param vin2 Second source rvec.
 * @param vout Destination rvec.
 */
void pbc_rvec_sub(const t_pbc *const pbc,
                  const rvec vin1,
                  const rvec vin2,
                  rvec vout) {
    if (pbc) {
        pbc_dx_aiuc(pbc, vin1, vin2, vout);
    }
    else {
        rvec_sub(vin1, vin2, vout);
    }
}
}   // namespace

bool Ilves::constr_dd() const {
#ifdef GMX_MPI
    return DOMAINDECOMP(cr) && cr->dd->constraints != nullptr;
#else
    return false;
#endif
}

void Ilves::omp_error_reduction(real &gerror, real perror) const {
    #pragma omp for reduction(max:gerror) nowait
    for (int t = 0; t < nthreads; ++t) {
        gerror = std::max(gerror, perror);
    }
}

void Ilves::dd_error_reduction(real &error) {
#ifdef GMX_MPI
    MPI_Allreduce(MPI_IN_PLACE, &error, 1, GMX_MPI_REAL, MPI_MAX, cr->dd->mpi_comm_all);
#endif
}

void Ilves::dd_comm_shared_xprime(const matrix box, const ArrayRef<RVec> xprime) {
#ifdef GMX_MPI
    dd_move_x_constraints(cr->dd, box, xprime, ArrayRef<RVec>(), FALSE);
#endif
}

real Ilves::make_rhs_scalar(const t_pbc *const pbc,
                            const ArrayRef<const RVec> x,
                            const ArrayRef<const RVec> xprime,
                            const bool compute_x_ab,
                            real *const rhs,
                            const int gstart,
                            const int gend,
                            const int lstart) {

    const bool xprime_ab_empty = xprime_ab.back().empty();

    real rel = 0;

    for (int grow = gstart, lrow = lstart; grow < gend; ++grow, ++lrow) {
        // Isolate the atoms which take part in bond row
        const int a = mol->bonds.atom1[grow];
        const int b = mol->bonds.atom2[grow];

        // Compute the vectors from atom a to atom b.
        if (compute_x_ab) {
            RVec rab;
            pbc_rvec_sub(pbc, x[b], x[a], rab);
            for (int d = 0; d < DIM; ++d) {
                x_ab[d][grow] = rab[d];
            }
        }
        RVec rcd;
        pbc_rvec_sub(pbc, xprime[b], xprime[a], rcd);
        if (!xprime_ab_empty) {
            for (int d = 0; d < DIM; ++d) {
                xprime_ab[d][grow] = rcd[d];
            }
        }

        // Compute the square of the length of r
        const real scalar = ::iprod(rcd, rcd);

        // Compute the constraint violation
        rhs[lrow] = 0.5 * (scalar - mol->bonds.sigma2[grow]);

        // Update the relative error
        rel = std::max(rel, std::abs(rhs[lrow]) * mol->bonds.invsigma2[grow]);
    }
    // Return the largest relative (square) bond length violation
    return rel;
}

#if GMX_SIMD_HAVE_REAL
real Ilves::make_rhs_simd(const real *const pbc_simd,
                          const ArrayRef<const RVec> x,
                          const ArrayRef<const RVec> xprime,
                          const bool compute_x_ab,
                          real *const g,
                          const int gstart,
                          const int gend) {

    const bool xprime_ab_empty = xprime_ab.back().empty();

    const real *x_data = reinterpret_cast<const real *>(as_rvec_array(x.data()));
    const real *
        xprime_data = reinterpret_cast<const real *>(as_rvec_array(xprime.data()));

    SimdReal vrel = 0;

    for (int grow = gstart, lrow = 0; grow < gend;
         grow += GMX_SIMD_REAL_WIDTH, lrow += GMX_SIMD_REAL_WIDTH) {

    #if GMX_DOUBLE
        auto va = simdLoadU(mol->bonds.atom1.data() + grow, SimdDInt32Tag());
        auto vb = simdLoadU(mol->bonds.atom2.data() + grow, SimdDInt32Tag());
    #else
        auto va = simdLoadU(mol->bonds.atom1.data() + grow, SimdFInt32Tag());
        auto vb = simdLoadU(mol->bonds.atom2.data() + grow, SimdFInt32Tag());
    #endif

        SimdReal va_x, va_y, va_z;
        SimdReal vb_x, vb_y, vb_z;
        SimdReal vrab_x, vrab_y, vrab_z;
        if (compute_x_ab) {
            // Compute the vectors from atom a to atom b.
            // pbc_rvec_sub(pbc, x[b], x[a], x_ab[grow]);
            gatherLoadUBySimdIntTranspose<DIM>(x_data, va, &va_x, &va_y, &va_z);
            gatherLoadUBySimdIntTranspose<DIM>(x_data, vb, &vb_x, &vb_y, &vb_z);

            vrab_x = vb_x - va_x;
            vrab_y = vb_y - va_y;
            vrab_z = vb_z - va_z;

            pbc_correct_dx_simd(&vrab_x, &vrab_y, &vrab_z, pbc_simd);

            storeU(x_ab[XX].data() + grow, vrab_x);
            storeU(x_ab[YY].data() + grow, vrab_y);
            storeU(x_ab[ZZ].data() + grow, vrab_z);
        }

        // pbc_rvec_sub(pbc, xprime[b], xprime[a], xprime_ab[grow]);
        gatherLoadUBySimdIntTranspose<DIM>(xprime_data, va, &va_x, &va_y, &va_z);
        gatherLoadUBySimdIntTranspose<DIM>(xprime_data, vb, &vb_x, &vb_y, &vb_z);

        vrab_x = vb_x - va_x;
        vrab_y = vb_y - va_y;
        vrab_z = vb_z - va_z;

        pbc_correct_dx_simd(&vrab_x, &vrab_y, &vrab_z, pbc_simd);

        if (!xprime_ab_empty) {
            storeU(xprime_ab[XX].data() + grow, vrab_x);
            storeU(xprime_ab[YY].data() + grow, vrab_y);
            storeU(xprime_ab[ZZ].data() + grow, vrab_z);
        }

        // Compute the square of the length of r
        // const auto scalar = iprod(xprime_ab[grow], xprime_ab[grow]);
        SimdReal vscalar = iprod(vrab_x, vrab_y, vrab_z, vrab_x, vrab_y, vrab_z);

        // Compute the constraint violation
        SimdReal vsigma2 = loadU<SimdReal>(mol->bonds.sigma2.data() + grow);
        SimdReal vinvsigma2 = loadU<SimdReal>(mol->bonds.invsigma2.data() + grow);

        vscalar = 0.5 * (vscalar - vsigma2);

        store(g + lrow, vscalar);

        vrel = max(vrel, abs(vscalar) * vinvsigma2);
    }

    alignas(GMX_SIMD_ALIGNMENT) real vrel_mem[GMX_SIMD_REAL_WIDTH];
    store(vrel_mem, vrel);

    real rel = 0;
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; ++i) {
        rel = std::max(rel, vrel_mem[i]);
    }

    return rel;
}
#endif

real Ilves::make_rhs(const int partition,
                     const t_pbc *const pbc,
                     const real *const pbc_simd,
                     const ArrayRef<const RVec> x,
                     const ArrayRef<const RVec> xprime,
                     const bool compute_x_ab) {

    auto &pdata = schur_solver->part_data[partition];

    // Nullify the schur entries of the rhs.
    #pragma omp simd
    for (int row = pdata.local_rows; row < pdata.local_rows + pdata.schur_rows;
         ++row) {
        pdata.rhs[row] = 0;
    }

    // Initialize the largest relative error
    real rel = 0;
    real new_rel = 0;

    int scalar_gstart = pdata.part[0];
    int scalar_lstart = 0;

#if GMX_SIMD_HAVE_REAL
    int simd_gstart = pdata.part[0];
    int simd_gend = pdata.part[0] +
                    ((pdata.part[1] - pdata.part[0]) / GMX_SIMD_REAL_WIDTH) *
                        GMX_SIMD_REAL_WIDTH;

    new_rel = make_rhs_simd(pbc_simd,
                            x,
                            xprime,
                            compute_x_ab,
                            pdata.rhs.data(),
                            simd_gstart,
                            simd_gend);
    rel = std::max(rel, new_rel);

    scalar_gstart = simd_gend;
    scalar_lstart = simd_gend - simd_gstart;
#endif
    new_rel = make_rhs_scalar(pbc,
                              x,
                              xprime,
                              compute_x_ab,
                              pdata.rhs.data(),
                              scalar_gstart,
                              pdata.part[1],
                              scalar_lstart);
    rel = std::max(rel, new_rel);

    // Return the largest relative (square) bond length violation
    return rel;
}

void Ilves::make_lhs_scalar(const int partition,
                            const std::array<std::vector<real, AlignedAllocator<real>>,
                                             DIM> &xab1,
                            const std::array<std::vector<real, AlignedAllocator<real>>,
                                             DIM> &xab2,
                            const int lrowstart) {

    auto &pdata = schur_solver->part_data[partition];
    const auto &weights = part_lhs_weights[partition];

    const int lrowend = pdata.lhs.size();
    for (int lrow = lrowstart; lrow < lrowend; ++lrow) {
        const int grow = pdata.grows[lrow];
        const int gcol = pdata.gcols[lrow];

        // Compute the inner product scalar.
        const real scalar = xab1[XX][grow] * xab2[XX][gcol] +
                            xab1[YY][grow] * xab2[YY][gcol] +
                            xab1[ZZ][grow] * xab2[ZZ][gcol];

        pdata.lhs[lrow] = weights[lrow] * scalar;
    }
}

#if GMX_SIMD_HAVE_REAL
void Ilves::make_lhs_simd(const int partition,
                          const std::array<std::vector<real, AlignedAllocator<real>>,
                                           DIM> &xab1,
                          const std::array<std::vector<real, AlignedAllocator<real>>,
                                           DIM> &xab2,
                          const int lrowend) {
    auto &pdata = schur_solver->part_data[partition];

    const auto *weights = part_lhs_weights[partition].data();
    auto *lhs = pdata.lhs.data();

    const auto *xab1_XX = xab1[XX].data();
    const auto *xab1_YY = xab1[YY].data();
    const auto *xab1_ZZ = xab1[ZZ].data();

    const auto *xab2_XX = xab2[XX].data();
    const auto *xab2_YY = xab2[YY].data();
    const auto *xab2_ZZ = xab2[ZZ].data();

    const auto *grows = pdata.grows.data();
    const auto *gcols = pdata.gcols.data();

    for (int lrow = 0; lrow < lrowend; lrow += GMX_SIMD_REAL_WIDTH) {
    #if GMX_DOUBLE
        auto vgrow = simdLoad(grows + lrow, SimdDInt32Tag());
        auto vgcol = simdLoad(gcols + lrow, SimdDInt32Tag());
    #else
        auto vgrow = simdLoad(grows + lrow, SimdFInt32Tag());
        auto vgcol = simdLoad(gcols + lrow, SimdFInt32Tag());
    #endif

        SimdReal va_x, va_y, va_z;
        gatherLoadUBySimdIntTranspose<1>(xab1_XX, vgrow, &va_x);
        gatherLoadUBySimdIntTranspose<1>(xab1_YY, vgrow, &va_y);
        gatherLoadUBySimdIntTranspose<1>(xab1_ZZ, vgrow, &va_z);
        SimdReal vb_x, vb_y, vb_z;
        gatherLoadUBySimdIntTranspose<1>(xab2_XX, vgcol, &vb_x);
        gatherLoadUBySimdIntTranspose<1>(xab2_YY, vgcol, &vb_y);
        gatherLoadUBySimdIntTranspose<1>(xab2_ZZ, vgcol, &vb_z);

        // Compute the inner product scalar.
        SimdReal vscalar = iprod(va_x, va_y, va_z, vb_x, vb_y, vb_z);

        SimdReal vweight = load<SimdReal>(weights + lrow);

        store(lhs + lrow, vweight * vscalar);
    }
}
#endif

void Ilves::make_lhs(const int partition,
                     const std::array<std::vector<real, AlignedAllocator<real>>, DIM> &xab1,
                     const std::array<std::vector<real, AlignedAllocator<real>>, DIM> &xab2) {

    auto &pdata = schur_solver->part_data[partition];

    int scalar_lrowstart = 0;

#if GMX_SIMD_HAVE_REAL
    const int simd_lrowend = (pdata.lhs.size() / GMX_SIMD_REAL_WIDTH) *
                             GMX_SIMD_REAL_WIDTH;

    make_lhs_simd(partition, xab1, xab2, simd_lrowend);

    scalar_lrowstart = simd_lrowend;
#endif
    make_lhs_scalar(partition, xab1, xab2, scalar_lrowstart);
}

void Ilves::update_current_lagr(const int partition, const bool first_time) {
    auto &pdata = schur_solver->part_data[partition];

    if (first_time) {
        std::swap(current_lagr[partition], pdata.rhs);
    }
    else {
        auto *rhs_data = pdata.rhs.data();
        auto *current_lagr_data = current_lagr[partition].data();

        #pragma omp simd aligned(rhs_data, current_lagr_data: GMX_SIMD_ALIGNMENT)
        for (int lrow = 0; lrow < pdata.local_rows; ++lrow) {
            current_lagr_data[lrow] += rhs_data[lrow];
        }
    }
}

void Ilves::update_positions(const int partition, const ArrayRef<RVec> xprime) const {
    const auto &pdata = schur_solver->part_data[partition];

    for (int grow = pdata.part[0], lrow = 0; grow < pdata.part[1]; ++grow, ++lrow) {
        const int a = mol->bonds.atom1[grow];
        const int b = mol->bonds.atom2[grow];

        // Update the components of y.
        const real rhs_a = pdata.rhs[lrow] * mol->atoms.invmass[a];
        const real rhs_b = pdata.rhs[lrow] * mol->atoms.invmass[b];
        for (int d = 0; d < DIM; ++d) {
            xprime[a][d] += rhs_a * x_ab[d][grow];
            xprime[b][d] -= rhs_b * x_ab[d][grow];
        }
    }
}

void Ilves::update_virial(const int partition, tensor virial) const {
    const auto &pdata = schur_solver->part_data[partition];

    for (int grow = pdata.part[0], lrow = 0; grow < pdata.part[1]; ++grow, ++lrow) {
        for (int d1 = 0; d1 < DIM; ++d1) {
            const real tmp = -current_lagr[partition][lrow] * x_ab[d1][grow];

            for (int d2 = 0; d2 < DIM; d2++) {
                virial[d1][d2] -= tmp * x_ab[d2][grow];
            }
        }
    }
}

void Ilves::update_velocities(const int partition,
                              const ArrayRef<RVec> vprime,
                              const real invdt) const {
    const auto &pdata = schur_solver->part_data[partition];

    for (int grow = pdata.part[0], lrow = 0; grow < pdata.part[1]; ++grow, ++lrow) {
        const int a = mol->bonds.atom1[grow];
        const int b = mol->bonds.atom2[grow];

        // Update the components of y.
        const real lgr_invdt = current_lagr[partition][lrow] * invdt;
        const real lgr_invdt_a = lgr_invdt * mol->atoms.invmass[a];
        const real lgr_invdt_b = lgr_invdt * mol->atoms.invmass[b];
        for (int d = 0; d < DIM; ++d) {
            vprime[a][d] += lgr_invdt_a * x_ab[d][grow];

            vprime[b][d] -= lgr_invdt_b * x_ab[d][grow];
        }
    }
}

bool Ilves::disjoint_mol(const int submol_max_size) const {
    auto &graph = mol->bonds.graph;

    std::vector<bool> visited(graph.num_nodes(), false);

    for (int bond = 0; bond < graph.num_nodes(); ++bond) {
        if (visited[bond]) {
            continue;
        }

        visited[bond] = true;

        int nneighs = 1;

        // Visit all the neighbors of bond.
        std::function<void(const int)> count_neighs = [&](const int node) -> void {
            for (int k = graph.xadj[node]; k < graph.xadj[node + 1]; ++k) {
                const int neigh = graph.adj[k];

                if (visited[neigh]) {
                    continue;
                }

                visited[neigh] = true;

                ++nneighs;

                if (nneighs > submol_max_size) {
                    return;
                }

                count_neighs(neigh);
            }
        };

        count_neighs(bond);

        if (nneighs > submol_max_size) {
            return false;
        }
    }

    return true;
}

void Ilves::make_weights() {
    /*
       Description of the construction of the weights:

       The following example explains the origins of the weights and how to
       compute them.

       EXAMPLE:

       A molecule with m = 6 atoms and n = 5 bonds. The atoms are numbered 0
       through 5, and the bonds are numbered 0 through 4. It is irrelevant if
       this molecule exists in the real world or not :)


        0                 4           The bonds are given by the array
         \               /
         (0)           (3)            bonds = {0,2; 1,2; 2,3; 3,4; 3,5}
           \           /
            2 ------- 3               Each bond involves 2 atoms. Bond 2 is
           /    (2)    \              between atoms 2 and 3.
         (1)           (4)
         /               \
       1                  5

       The adjacency graph for the bonds is the graph


          0   3              matrix A is 5 by 5
         / \ / \
        |   2   |            |xxx  |
         \ / \ /             |xxx  |
          1   4              |xxxxx|
                             |  xxx|
                             |  xxx|

       This is how the graph is encoded

        adj  = {0,1,2; 0,1,2; 0,1,2,3,4; 2,3,4; 2,3,4}     (17 entries)
       xadj  = {0, 3, 6, 11, 14, 17}                       ( 6 entries)

       Please note that there are n+1 entries in XADJ, and the last entry point
       just beyond the end of adj. Therefore XADJ[n+1] is the number of nonzero
       entries in the matrix.

       NOTATION:

         1) We write rij or r(i,j) for the vector from atom i to j.
         2) We write <x,y> for the scalar product between the vectors x and y.

       In the above example bond number 2 involves atoms 2 and 3 and is a
       (mathematical) constraint of the type

                  0.5*(||r23||^2 - (bond length)^2) = 0

       The factor 0.5 is only included to give the Jacobian of the constraint
       function a form which I (CCKM) consider aesthetically pleasing.

       Now, our matrix of the form

                  A =  Dg(r)*inv(M)*Dg(s)'

       where r and s are vectors describing two different configurations of
       the m atoms, so r and s each have 3m components each.

       Below is a table of the nonzero entries of the matrix A for our current
       example:

                   bond = {0,2,1,2,2,3,3,4,3,5}

       i    j      entry A(i,j)
       ---------------------------------------------------------------
       0    0     (invmass(0)+invmass(2))*<r02,s03>
       1    0     +invmass(2)*<r12,s02>
       2    0     -invmass(2)*<r23,s02>

       0    1     +invmass(2)*<r02,s12>
       1    1     (invmass(1)+invmass(2))*<r12,s12>
       2    1     -invmass(2)*<r23,s12>

       0    2     -invmass(2)*<r02,s23>
       1    2     -invmass(2)*<r12,s23>
       2    2     (invmass(2)+invmass(3))*<r23,s23>
       3    2     -invmass(3)*<r34,s23>
       4    2     -invmass(3)*<r35,s23>

       2    3     -invmass(3)*<r23,s34>
       3    3     (invmass(3)+invmass(4))*<r34,s34>
       4    3     +invmass(3)*<r35,s34>

       2    4     -invmass(3)*<r23,s35>
       3    4     +invmass(3)*<r34,s35>
       4    4     (+invmass(3)+invmass(5))*<r35,s35>

       This table is very carefully constructed! Please note the following:

       a) Reading the table from the top to the bottom takes us through the
       nonzero entries of A in column major format, exactly as matrices are
       stored in LAPACK. Moreover, we are writing to main memory in an order
       which is strictly increasing.

       b) For each column of the matrix A we deal with a FIXED vector s(a,b),
       rather than both s(a,b) and s(b,a) = -s(a,b).

       c) The order of the indices as in r(a,b) for a bond k, is EXACTLY the
       order in which the atoms which partake in bond k are given in the bond
       table.

       EXAMPLE: As noted above bond 2 is a bond between atoms 2 and 3. The bond
       table lists atom 2 BEFORE atom 3. Therefore, we use the vector r23,
       rather than the (equvalent) vector r32 = -r23

       d) The diagonal entries are different from the off diagonal entries
       because the weights are different.

       In order to efficiently generate the matrix A we precompute the following
       weights

       i    j     weight
       ----------------------------------------------
       0    0     +invmass(0)+invmass(2)
       1    0     +invmass(2)
       2    0     -invmass(2)

       0    1     +invmass(2)
       1    1     +invmass(1)+invmass(2)
       2    1     -invmass(2)

       0    2     -invmass(2)
       1    2     -invmass(2)
       2    2     +invmass(2)+invmass(3)
       3    2     -invmass(3)
       4    2     -invmass(3)

       2    3     -invmass(3)
       3    3     +invmass(3)+invmass(4)
       4    3     +invmass(3)

       2    4     -invmass(3)
       3    4     +invmass(3)
       4    4     +invmass(3)+invmass(5)

       In short, the weights includes all the relevant information about the
       signs as well as the masses of the atoms.

       Given vectors r and s with a 3m components as well as the bond graph, we
       can now generate the matrix A one entry at a time, moving through RAM
       memory in a strictly increasing order.

    */
    const int nparts = schur_solver->part_data.size();
    part_lhs_weights.resize(nparts);

    for (int p = 0; p < nparts; ++p) {
        const auto &pdata = schur_solver->part_data[p];
        const auto &pmatrix = pdata.fill_matrix;

        auto &pweights = part_lhs_weights[p];
        pweights.resize(pmatrix.num_edges());

        for (int i = 0; i < pmatrix.num_edges(); ++i) {
            const int row = pdata.grows[i];
            const int col = pdata.gcols[i];

            // Isolate the atoms which partake in the rowth bond
            const int arow1 = mol->bonds.atom1[row];
            const int arow2 = mol->bonds.atom2[row];

            if (pdata.is_fillin[i]) {
                pweights[i] = 0;
            }
            else if (row != col) {
                // Isolate the atoms which partake in the colth bond
                const int acol1 = mol->bonds.atom1[col];
                const int acol2 = mol->bonds.atom2[col];

                // Determine the atom which is common to bond row and
                // bond col
                const int common = ((arow1 == acol1) || (arow1 == acol2)) ? arow1
                                                                          : arow2;

                pweights[i] = mol->atoms.invmass[common];

                // Determine the appropriate sign
                if ((arow1 == acol2) || (arow2 == acol1)) {
                    /* You should reverse the order the atoms for one of
                    the two bonds, but this impractical, so we just
                    reverse the sign of the weight
                    */
                    pweights[i] = -pweights[i];
                }
            }
            else {
                /* This is a diagonal entry.
                    The weight is special, but the order of the atoms in the
                    bond list is irrelevant. Yes, you could flip the order of
                    one pair of atoms, but then you would be compelled to
                    flip the order of the second pair, and so you would
                    change sign twice.
                */
                pweights[i] = mol->atoms.invmass[arow1] + mol->atoms.invmass[arow2];
            }
        }
    }
}

}   // namespace gmx
