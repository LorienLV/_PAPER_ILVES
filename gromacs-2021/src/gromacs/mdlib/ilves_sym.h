#pragma once

#include "gromacs/utility/real.h"
#include "ilves.h"

/**
 * ILVES implementation that uses the symmetric matrix approach.
 *
 */

namespace gmx {

class IlvesSym : public Ilves {
public:
    /**
     * Initializes the ILVES solver for constraining the bonds described by
     * idef.
     *
     * @param cr MPI communication structure.
     * @param idef The interaction definitions provided by GROMACS.
     * @param invmass A pointer to the array of inverse mass of each atom.
     * ILVES keeps a pointer to this array, so it must be kept alive.
     * @param nthreads Number of threads to use.
     * @param durations Reference to detailed timing information.
     */
    IlvesSym(const t_commrec *cr,
             const InteractionDefinitions &idef,
             const real *invmass,
             int threads,
             std::vector<std::map<std::string, std::chrono::duration<double, std::milli>>>
                 &durations);

    /**
     * Try to solve the bond constraints of MOL in at most MAXIT iterations of
     * Newton's method with a tolerance for each atom of TOL.
     *
     * @param x Positions of the atoms before computing the external forces,
     * with the following format:
     *
     * atom 1     atom 2     atom 3
     * x, y, z,   x, y, z,   x, y, z
     *
     * @param xprime Positions of the atoms after computing the external forces,
     * with the same format as x. When returning, it will contain the final
     * position of each atom.
     * @param v The atoms velocities. When returning, it will contain the final
     * velocity of each atom.
     * @param tol Maximum error allowed in each atom position.
     * @param deltat The time-step in picoseconds.
     * @param constraint_virial Update the virial?
     * @param virial sum r x m delta_r (virial)
     * @param constraint_velocities Update the velocities?
     * @param pbc The PBC (container) information. Null if there is no PBC
     * information.
     * @param compute_fep Compute the free energy?
     * @param fep_lambda FEP lambda.
     * @param fep_force The FEP force.
     * @param box Used for pbc in MPI runs.
     * @param perf_stats Performance statistics.
     * @return An std::pair. The first element is true if the constraints has
     * been satisfied, false otherwise. The second element is the number of
     * iterations performed by the solver.
     */
    std::pair<bool, int> solve(ArrayRef<const RVec> x,
                               ArrayRef<RVec> xprime,
                               ArrayRef<RVec> vprime,
                               real tol,
                               real deltat,
                               bool constraint_virial,
                               tensor virial,
                               bool constraint_velocities,
                               const t_pbc *pbc,
                               bool compute_fep,
                               real fep_lambda,
                               real &fep_force,
                               const matrix box,
                               t_nrnb *perf_stats) override;
};
}   // namespace gmx