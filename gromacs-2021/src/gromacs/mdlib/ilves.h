#pragma once

#include <chrono>
#include <map>

#include "gromacs/domdec/domdec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"
#include "gromacs/simd/simd.h"

#include "molecule.h"
#include "schur_linear_solver.h"

#define PRINT_FUNCTION_DURATIONS 0

/**
 * ILVES abstract class. This class constains all the methods required by a
 * child implementation to solve the constraints of a molecule.
 *
 */

namespace gmx {

class Ilves {
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
     * @param upper_tri ILVES will use a symmetric matrix to solve the
     * constraints. Otherwise, ILVES will use a structurally symmetric matrix.
     * @param durations Reference to detailed timing information.
     */
    Ilves(const t_commrec *cr,
          const InteractionDefinitions &idef,
          const real *invmass,
          int nthreads,
          bool upper_tri,
          std::vector<std::map<std::string, std::chrono::duration<double, std::milli>>>
              &durations);

    /**
     * Destroy ILVES. This is an abstract class, so we define it virtual to
     * avoid problems.
     *
     */
    virtual ~Ilves() = default;

    /**
     * This is a virtual function.
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
    virtual std::pair<bool, int> solve(ArrayRef<const RVec> x,
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
                                       t_nrnb *perf_stats) = 0;

protected:
    int nthreads;   // Number of threads to use.

    static constexpr int maxit = 16;   // Maximum allowed iterations.

    std::unique_ptr<Molecule> mol;   // The molecule structure.

    // Parallel linear solver.
    std::unique_ptr<SchurLinearSolver> schur_solver;
    // Weights of the entries of the lhs.
    std::vector<std::vector<real, AlignedAllocator<real>>> part_lhs_weights;

    // The current approximation of the Lagrange multipliers. One vector per
    // partition.
    std::vector<std::vector<real, AlignedAllocator<real>>> current_lagr;

    // x_ab[XX/YY/ZZ][i] contains the vector from atom a to b (using x),
    // the atoms that are part of the ith bond.
    std::array<std::vector<real, AlignedAllocator<real>>, DIM> x_ab;
    // x_ab[XX/YY/ZZ][i] contains the vector from atom a to b (using xprime),
    // the atoms that are part of the ith bond.
    std::array<std::vector<real, AlignedAllocator<real>>, DIM> xprime_ab;

    // Counting execution time.
    // Fine grained timing information.
    std::vector<std::map<std::string, std::chrono::duration<double, std::milli>>> &durations;

    /**
     * True if this simulation uses more than one MPI rank for constraints.
     *
     * @return True if this simulation uses more than one MPI rank for
     * constraints.
     */
    bool constr_dd() const;

    /**
     * OpenMP max reduction of private taus into the global tau. We need this
     * function since some compilers produce slower code when performing a
     * reduction into a pointer variable.
     * See https://stackoverflow.com/questions/76480632
     *
     * @param gerror The global shared error.
     * @param perror The private error.
     */
    void omp_error_reduction(real &gerror, real perror) const;

    /**
     * MPI max reduction of error between all the ranks of the domain
     * decomposition.
     *
     * @param error The error of the rank.
     */
    void dd_error_reduction(real &error);

    /**
     * Wrapper for dd_move_x_constraints (communication of non-local atoms).
     *
     * @param box Used for pbc.
     * @param xprime Positions of the atoms.
     */
    void dd_comm_shared_xprime(const matrix box, ArrayRef<RVec> xprime);

    /**
     * Computes the part of the right-hand side of the linear system between
     * GSTART and GEND - 1 rows. Returns the largest relative bond length
     * violation of that part.
     *
     * @param pbc The PBC (container) information. Null if there is no PBC.
     * @param x Positions of the atoms before computing the external forces.
     * @param xprime Updated positions of the atoms.
     * @param compute_x_ab Update x_ab?
     * @param rhs Pointer to the local right hand side of the linear system.
     * @param gstart First global row id of the linear system.
     * @param gend Last global row id + 1 of the linear system.
     * @param lstart First local row id of the linear system.
     * @return The largest relative bond length violation of the part.
     */
    real make_rhs_scalar(const t_pbc *pbc,
                         ArrayRef<const RVec> x,
                         ArrayRef<const RVec> xprime,
                         bool compute_x_ab,
                         real *rhs,
                         int gstart,
                         int gend,
                         int lstart);

    /**
     * Computes the part of the right-hand side of the linear system between
     * GSTART and GEND - 1 rows using SIMD intrinsics. Returns the largest
     * relative bond length violation of the part. GEND - GSTART must be
     * a multiple of GMX_SIMD_REAL_WIDTH.
     *
     * @param pbc The PBC (container) information. Null if there is no PBC.
     * @param x Positions of the atoms before computing the external forces.
     * @param xprime Updated positions of the atoms.
     * @param compute_x_ab Update x_ab?
     * @param rhs Pointer to the local right hand side of the linear system.
     * @param gstart First global row id of the linear system.
     * @param gend Last global row id + 1 of the linear system.
     * @return The largest relative bond length violation of the part.
     */
#if GMX_SIMD_HAVE_REAL
    real make_rhs_simd(const real *pbc_simd,
                       ArrayRef<const RVec> x,
                       ArrayRef<const RVec> xprime,
                       bool compute_x_ab,
                       real *rhs,
                       int gstart,
                       int gend);
#endif

    /**
     * Computes the right-hand side of the partition PARTITION linear system.
     * Returns the largest relative bond length violation of the partition.
     *
     * @param partition Partition id.
     * @param pbc The PBC (container) information. Null if there is no PBC.
     * @param pbc_simd The PBC (container) information for SIMD execution.
     * Null if there is no PBC.
     * @param x Positions of the atoms before computing the external forces.
     * @param xprime Updated positions of the atoms.
     * @param compute_x_ab Update x_ab?
     * @return The largest relative bond length violation of the partition.
     */
    real make_rhs(int partition,
                  const t_pbc *pbc,
                  const real *pbc_simd,
                  ArrayRef<const RVec> x,
                  ArrayRef<const RVec> xprime,
                  bool compute_x_ab);

    /**
     * Computes the left-hand side of the partition PARTITION linear system
     * starting at the lrowstarth local row.
     *
     * @param partition Partition id.
     * @param xab1 An array of vectors containing the vector from atom a to b
     * (xprime_ab or x_ab).
     * @param xab2 An array of vectors containing the vector from atom a to b
     * (xprime_ab or x_ab).
     */
    void
    make_lhs_scalar(int partition,
                    const std::array<std::vector<real, AlignedAllocator<real>>, DIM> &xab1,
                    const std::array<std::vector<real, AlignedAllocator<real>>, DIM> &xab2,
                    int lrowstart);

    /**
     * Computes the left-hand side of the partition PARTITION linear system
     * between the 0th and the lrowendth - 1 local rows using SIMD intrinsics.
     *
     * @param partition Partition id.
     * @param xab1 An array of vectors containing the vector from atom a to b
     * (xprime_ab or x_ab).
     * @param xab2 An array of vectors containing the vector from atom a to b
     * (xprime_ab or x_ab).
     */
#if GMX_SIMD_HAVE_REAL
    void
    make_lhs_simd(int partition,
                  const std::array<std::vector<real, AlignedAllocator<real>>, DIM> &xab1,
                  const std::array<std::vector<real, AlignedAllocator<real>>, DIM> &xab2,
                  int lrowend);
#endif

    /**
     * Computes the left-hand side of the partition PARTITION linear system.
     *
     * @param partition Partition id.
     * @param xab1 An array of vectors containing the vector from atom a to b
     * (xprime_ab or x_ab).
     * @param xab2 An array of vectors containing the vector from atom a to b
     * (xprime_ab or x_ab).
     */
    void make_lhs(int partition,
                  const std::array<std::vector<real, AlignedAllocator<real>>, DIM> &xab1,
                  const std::array<std::vector<real, AlignedAllocator<real>>, DIM> &xab2);

    /**
     * Updates the vector of current lagrange multipliers of the partition
     * PARTITION with the solution of the PARTITION linear system.
     *
     * @param partition Partition id.
     * @param first_time Is this function called in the first iteration of
     * ILVES? (this information is used for performance optimization).
     */
    void update_current_lagr(int partition, bool first_time);

    /**
     * For each bond b, that joins atoms a and b:
     *     XPRIME[a] += lagr_correction[b] * (X>[b] - X>[a]) * (1 /
     * mass[a]); XPRIME[b] -= lagr_correction[b] * (X>[b] - X>[a]) * (1 /
     * mass[a]);
     *
     * The function only updates the positions atoms connected by bonds of the
     * partition PARTITION.
     *
     * @param partition Partition id.
     * @param xprime Positions of the atoms after computing the external
     * forces (current positions of the atoms).
     */
    void update_positions(int partition, ArrayRef<RVec> xprime) const;

    /**
     * For each bond b that joins atoms a and b:
     *
     * for d1 in DIM
     *     for d2 in DIM
     *         VIRIAL[d1][d2] -= -lagr[b] * (X>[b] - X>[a])[d1] * (X>[b] -
     * X>[a])[d2]
     *
     * The function only updates the virial with lagrange multipliers of
     * partition PARTITION.
     *
     * @param partition Partition id.
     * @param virial The virial (a DIMxDIM matrix).
     */
    void update_virial(int partition, tensor virial) const;

    /**
     * For each bond b that joins atoms a and b:
     *     V[a] += lagr[b] * (X>[b] - X>[a]) * (1 / mass[a]) * INVDT;
     *     V[b] -= lagr[b] * (X>[b] - X>[a]) * (1 / mass[a]) * INVDT;
     *
     * The function only updates the velocities of the atoms connected by bonds
     * of the partition PARTITION.
     * @param partition Partition id.
     * @param vprime The atoms velocities.
     * @param invdt Inverse of the time-step in picoseconds.
     */
    void update_velocities(int partition, ArrayRef<RVec> vprime, real invdt) const;

private:
    // Domain decomposition (MPI only).
    const t_commrec *const cr;
    // Range of atoms received from other MPI ranks.
    // std::pair<int, int> non_local_atoms_range;

    /**
     * Return true if the molecule contains multiple submolecules with at most
     * submol_max_size atoms.
     *
     * @return bool True if the molecule contains multiple submolecules with at
     * most submol_max_size atoms. False otherwise.
     */
    bool disjoint_mol(int submol_max_size) const;

    /**
     * Compute the weight of each entry of the left-hand-size of the linear
     * system.
     *
     */
    void make_weights();
};
}   // namespace gmx