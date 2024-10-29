#include "ilves_asym.h"

#include <omp.h>

#include <cmath>
#include <utility>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_simd.h"
#include "gromacs/simd/simd.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"
#include "ilves.h"

namespace gmx {

IlvesAsym::IlvesAsym(const t_commrec *const cr,
                     const InteractionDefinitions &idef,
                     const real *const invmass,
                     const int nthreads,
                     std::vector<std::map<std::string, std::chrono::duration<double, std::milli>>>
                         &durations)
    : Ilves(cr, idef, invmass, nthreads, false, durations) {

    for (int d = 0; d < DIM; ++d) {
        xprime_ab[d].resize(mol->bonds.num);
    }
}

std::pair<bool, int> IlvesAsym::solve(const ArrayRef<const RVec> x,
                                      const ArrayRef<RVec> xprime,
                                      const ArrayRef<RVec> vprime,
                                      const real tol,
                                      const real deltat,
                                      const bool constraint_virial,
                                      tensor virial,
                                      const bool constraint_velocities,
                                      const t_pbc *const pbc,
                                      const bool compute_fep,
                                      const real fep_lambda,
                                      real &fep_force,
                                      const matrix box,
                                      t_nrnb *const perf_stats) {

    int numit = 0;

    // It might happen that all the partitions are disconnected, for example
    // when constraining h-bonds. In this case we can reduce the number of
    // synchronization points.
    const auto &schur_part = schur_solver->part_data.back().part;
    bool dep_par = (schur_part[1] - schur_part[0]) > 0;

    real taus[maxit + 1] = {0};

    #pragma omp parallel num_threads(nthreads)
    {
        const int tid = omp_get_thread_num();

        real *pbc_simd = nullptr;
#if GMX_SIMD_HAVE_REAL
        /* Convert the pbc struct for SIMD */
        alignas(GMX_SIMD_ALIGNMENT) real pbc_simd_data[9 * GMX_SIMD_REAL_WIDTH];
        pbc_simd = pbc_simd_data;

        set_pbc_simd(pbc, pbc_simd);
#endif

        /*
         * Compute g(x). Update x_ab and xprime_ab.
         */
#if PRINT_FUNCTION_DURATIONS
        auto start = std::chrono::system_clock::now();
#endif

        // General block g(x).
        real ptau = make_rhs(tid, pbc, pbc_simd, x, xprime, true);

        if (dep_par) {
            // Schur block g(x).
            #pragma omp master
            {
                real schur_tau = make_rhs(nthreads, pbc, pbc_simd, x, xprime, true);
                ptau = std::max(ptau, schur_tau);
            }
        }
#if PRINT_FUNCTION_DURATIONS
        auto end = std::chrono::system_clock::now();
        durations[tid]["make_rhs"] += end - start;
#endif

        omp_error_reduction(taus[0], ptau);
        #pragma omp barrier

#if PRINT_FUNCTION_DURATIONS
        start = std::chrono::system_clock::now();
#endif
        if (constr_dd()) {
            #pragma omp master
            {
                // MPI reduction of tau
                dd_error_reduction(taus[0]);
            }
            #pragma omp barrier
        }

        ptau = taus[0];

#if PRINT_FUNCTION_DURATIONS
        end = std::chrono::system_clock::now();
        durations[tid]["comm_mpi"] += end - start;
#endif

        // Do at most MAXIT Newton steps.
        for (int i = 0; i < maxit && std::isfinite(ptau) && tol < ptau; ++i) {
            #pragma omp master
            { ++numit; }

            /*
             * Compute A.
             */

#if PRINT_FUNCTION_DURATIONS
            start = std::chrono::system_clock::now();
#endif
            // General block A.
            make_lhs(tid, xprime_ab, x_ab);

            // Schur block A.
            if (dep_par) {
                #pragma omp master
                { make_lhs(nthreads, xprime_ab, x_ab); }
                // We need to wait for the shared lhs before starting the
                // factorization.
                #pragma omp barrier
            }
#if PRINT_FUNCTION_DURATIONS
            end = std::chrono::system_clock::now();
            durations[tid]["make_lhs"] += end - start;
#endif

            // No need to wait. The scheduling will be the same in the next
            // function.

            /*
             * Parallel LU factorization
             */
#if PRINT_FUNCTION_DURATIONS
            start = std::chrono::system_clock::now();
#endif
            schur_solver->LU_factor();
#if PRINT_FUNCTION_DURATIONS
            end = std::chrono::system_clock::now();
            durations[tid]["LU"] += end - start;
#endif

            // No need to wait. The scheduling will be the same in the next
            // function.

            /*
             * Parallel forward and backward substitution.
             */
#if PRINT_FUNCTION_DURATIONS
            start = std::chrono::system_clock::now();
#endif
            schur_solver->LU_solve();
#if PRINT_FUNCTION_DURATIONS
            end = std::chrono::system_clock::now();
            durations[tid]["forward + backward"] += end - start;
#endif

            // No need to wait. The scheduling will be the same in the next
            // function.

            /*
             * Update the positions of the atoms.
             */
#if PRINT_FUNCTION_DURATIONS
            start = std::chrono::system_clock::now();
#endif
            update_positions(tid, xprime);
#if PRINT_FUNCTION_DURATIONS
            end = std::chrono::system_clock::now();
            durations[tid]["update_positions"] += end - start;
#endif

            if (dep_par) {
                // Avoids that two threads modify the same atom position at the
                // same time.
                #pragma omp barrier

#if PRINT_FUNCTION_DURATIONS
                start = std::chrono::system_clock::now();
#endif

                #pragma omp master
                { update_positions(nthreads, xprime); }

#if PRINT_FUNCTION_DURATIONS
                end = std::chrono::system_clock::now();
                durations[tid]["update_positions"] += end - start;
#endif
                // Wait for updated positions.
                #pragma omp barrier
            }

#if PRINT_FUNCTION_DURATIONS
            start = std::chrono::system_clock::now();
#endif
            if (constr_dd()) {
                if (!dep_par) {
                    // We did not wait for all positions before.
                    #pragma omp barrier
                }

                // Communicate the corrected non-local coordinates (MPI).
                #pragma omp master
                { dd_comm_shared_xprime(box, xprime); }
                #pragma omp barrier
            }
#if PRINT_FUNCTION_DURATIONS
            end = std::chrono::system_clock::now();
            durations[tid]["comm_mpi"] += end - start;
#endif

            /*
             * Update the current approximation of the lagrange multipliers.
             */
#if PRINT_FUNCTION_DURATIONS
            start = std::chrono::system_clock::now();
#endif
            update_current_lagr(tid, i == 0);

            if (dep_par) {
                #pragma omp master
                { update_current_lagr(nthreads, i == 0); }
            }
#if PRINT_FUNCTION_DURATIONS
            end = std::chrono::system_clock::now();
            durations[tid]["update_lagr"] += end - start;
#endif

            // No need to wait. The scheduling will be the same in the next
            // function.

            /*
             * Compute g(x). Update xprime_ab only.
             */

            // General block g(x).
#if PRINT_FUNCTION_DURATIONS
            start = std::chrono::system_clock::now();
#endif
            ptau = make_rhs(tid, pbc, pbc_simd, x, xprime, false);

            if (dep_par) {
                // Schur block g(x).
                #pragma omp master
                {
                    real schur_tau = make_rhs(nthreads, pbc, pbc_simd, x, xprime, false);
                    ptau = std::max(ptau, schur_tau);
                }
            }
#if PRINT_FUNCTION_DURATIONS
            end = std::chrono::system_clock::now();
            durations[tid]["make_rhs"] += end - start;
#endif

            omp_error_reduction(taus[i + 1], ptau);
            #pragma omp barrier

#if PRINT_FUNCTION_DURATIONS
            start = std::chrono::system_clock::now();
#endif
            if (constr_dd()) {
                #pragma omp master
                {
                    // MPI reduction of tau
                    dd_error_reduction(taus[i + 1]);
                }
                #pragma omp barrier
            }

            ptau = taus[i + 1];

#if PRINT_FUNCTION_DURATIONS
            end = std::chrono::system_clock::now();
            durations[tid]["comm_mpi"] += end - start;
#endif
        }

        /*
         * Update the virial.
         */
#if PRINT_FUNCTION_DURATIONS
        start = std::chrono::system_clock::now();
#endif
        if (constraint_virial && numit > 0) {
            tensor pvirial;

            // Init the private virial to 0;
            for (int d1 = 0; d1 < DIM; ++d1) {
                for (int d2 = 0; d2 < DIM; ++d2) {
                    pvirial[d1][d2] = 0;
                }
            }

            update_virial(tid, pvirial);

            if (dep_par) {
                #pragma omp master
                { update_virial(nthreads, pvirial); }
            }

            // Update the global virial.
            for (int d1 = 0; d1 < DIM; ++d1) {
                for (int d2 = 0; d2 < DIM; ++d2) {
                    #pragma omp atomic
                    virial[d1][d2] += pvirial[d1][d2];
                }
            }
        }
#if PRINT_FUNCTION_DURATIONS
        end = std::chrono::system_clock::now();
        durations[tid]["update_virial"] += end - start;
#endif

        /*
         * Update the velocities.
         */
        if (constraint_velocities && numit > 0) {
#if PRINT_FUNCTION_DURATIONS
            start = std::chrono::system_clock::now();
#endif
            update_velocities(tid, vprime, 1 / deltat);
#if PRINT_FUNCTION_DURATIONS
            end = std::chrono::system_clock::now();
            durations[tid]["update_velocities"] += end - start;
#endif

            if (dep_par) {
                // Avoids that two threads modify the same atom velocity at the
                // same time.
                #pragma omp barrier

#if PRINT_FUNCTION_DURATIONS
                start = std::chrono::system_clock::now();
#endif

                #pragma omp master
                { update_velocities(nthreads, vprime, 1 / deltat); }

#if PRINT_FUNCTION_DURATIONS
                end = std::chrono::system_clock::now();
                durations[tid]["update_velocities"] += end - start;
#endif
            }
        }
    }

    /*
     * Update GROMACS performance statistics.
     */

    // The number of iterations.
    // TODO: Create eNR_ILVES.
    inc_nrnb(perf_stats, eNR_SHAKE, numit);
    // The number of constraints.
    // TODO: Create eNR_ILVES.
    inc_nrnb(perf_stats, eNR_SHAKE_RIJ, mol->bonds.num);
    if (constraint_velocities) {
        // The number of constrained velocities.
        inc_nrnb(perf_stats, eNR_CONSTR_V, mol->bonds.num * 2);
    }
    if (constraint_virial) {
        // The number of constrained virials.
        inc_nrnb(perf_stats, eNR_CONSTR_VIR, mol->bonds.num);
    }

    const real last_tau = taus[numit];
    return std::make_pair(std::isfinite(last_tau) && last_tau < tol, numit);
}

}   // namespace gmx
