/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief Defines the high-level constraint code.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "constr.h"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdlib/settle.h"
#include "gromacs/mdlib/shake.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/txtdump.h"

// ILVES
#include <chrono>
#include <iostream>
#include <iomanip>

#include "gromacs/simd/simd.h"

#include "omp.h"

#include "ilves.h"
#include "ilves_asym.h"
#include "ilves_sym.h"
// END ILVES


namespace gmx
{

/* \brief Impl class for Constraints
 *
 * \todo Members like md, idef are valid only for the lifetime of a
 * domain, which would be good to make clearer in the structure of the
 * code. It should not be possible to call apply() if setConstraints()
 * has not been called. For example, this could be achieved if
 * setConstraints returned a valid object with such a method.  That
 * still requires that we manage the lifetime of that object
 * correctly, however. */
class Constraints::Impl
{
public:
    Impl(const gmx_mtop_t&     mtop_p,
         const t_inputrec&     ir_p,
         pull_t*               pull_work,
         FILE*                 log_p,
         const t_commrec*      cr_p,
         const gmx_multisim_t* ms,
         t_nrnb*               nrnb,
         gmx_wallcycle*        wcycle_p,
         bool                  pbcHandlingRequired,
         int                   numConstraints,
         int                   numSettles);
    ~Impl();
    void setConstraints(gmx_localtop_t* top,
                        int             numAtoms,
                        int             numHomeAtoms,
                        real*           masses,
                        real*           inverseMasses,
                        bool            hasMassPerturbedAtoms,
                        real            lambda,
                        unsigned short* cFREEZE);
    bool apply(bool                      bLog,
               bool                      bEner,
               int64_t                   step,
               int                       delta_step,
               real                      step_scaling,
               ArrayRefWithPadding<RVec> x,
               ArrayRefWithPadding<RVec> xprime,
               ArrayRef<RVec>            min_proj,
               const matrix              box,
               real                      lambda,
               real*                     dvdlambda,
               ArrayRefWithPadding<RVec> v,
               bool                      computeVirial,
               tensor                    constraintsVirial,
               ConstraintVariable        econq);

    // ILVES
    std::pair<real, real> compute_errors(const t_pbc *const pbc,
                                            ArrayRef<RVec> xprime);
    void print_constr_detailed_durations();
    void mpi_start_printing();
    void mpi_finish_printing();
    // END ILVES

    //! The total number of constraints.
    int ncon_tot = 0;
    //! The number of flexible constraints.
    int nflexcon = 0;
    //! A list of atoms to constraints for each moleculetype.
    std::vector<ListOfLists<int>> at2con_mt;
    //! A list of atoms to settles for each moleculetype
    std::vector<std::vector<int>> at2settle_mt;
    //! LINCS data.
    Lincs* lincsd = nullptr; // TODO this should become a unique_ptr
    //! SHAKE data.
    std::unique_ptr<shakedata> shaked;
    //! SETTLE data.
    std::unique_ptr<SettleData> settled;
    //! The maximum number of warnings.
    int maxwarn = 0;
    //! The number of warnings for LINCS.
    int warncount_lincs = 0;
    //! The number of warnings for SETTLE.
    int warncount_settle = 0;
    //! The essential dynamics data.
    gmx_edsam* ed = nullptr;

    //! Thread-local virial contribution.
    tensor* threadConstraintsVirial = { nullptr };
    //! Did a settle error occur?
    bool* bSettleErrorHasOccurred = nullptr;

    //! Pointer to the global topology - only used for printing warnings.
    const gmx_mtop_t& mtop;
    //! Parameters for the interactions in this domain.
    const InteractionDefinitions* idef = nullptr;
    //! Total number of atoms.
    int numAtoms_ = 0;
    //! Number of local atoms.
    int numHomeAtoms_ = 0;
    //! Atoms masses.
    real* masses_;
    //! Inverse masses.
    real* inverseMasses_;
    //! If there are atoms with perturbed mass (for FEP).
    bool hasMassPerturbedAtoms_;
    //! FEP lambda value.
    real lambda_;
    //! Freeze groups data
    unsigned short* cFREEZE_;
    //! Whether we need to do pbc for handling bonds.
    bool pbcHandlingRequired_ = false;

    //! Logging support.
    FILE* log = nullptr;
    //! Communication support.
    const t_commrec* cr = nullptr;
    //! Multi-sim support.
    const gmx_multisim_t* ms = nullptr;
    //! Pulling code object, if any.
    pull_t* pull_work = nullptr;
    /*!\brief Input options.
     *
     * \todo Replace with IMdpOptions */
    const t_inputrec& ir;
    //! Flop counting support.
    t_nrnb* nrnb = nullptr;
    //! Tracks wallcycle usage.
    gmx_wallcycle* wcycle;

    // ILVES

    // 0 no
    // 1 yes
    // 2 ultra verbose
    int CONSTRAINTS_VERBOSE;

    // Number of steps that elapse between printing the constraints verbose
    // output
    int NSCONSTRAINTS_OUTPUT;

    std::unique_ptr<Ilves> ilves;

    int constr_total_calls = 0;
    int constr_total_iters = 0;
    std::chrono::duration<double> constr_duration = std::chrono::duration<double>(0.0);
    std::chrono::duration<double> constr_init_duration = std::chrono::duration<double>(0.0);
    std::vector<std::map<std::string,
                std::chrono::duration<double, std::milli>>> constr_detailed_durations;

    // END ILVES
};

Constraints::~Constraints() = default;

int Constraints::numFlexibleConstraints() const
{
    return impl_->nflexcon;
}

bool Constraints::havePerturbedConstraints() const
{
    const gmx_ffparams_t& ffparams = impl_->mtop.ffparams;

    for (size_t i = 0; i < ffparams.functype.size(); i++)
    {
        if ((ffparams.functype[i] == F_CONSTR || ffparams.functype[i] == F_CONSTRNC)
            && ffparams.iparams[i].constr.dA != ffparams.iparams[i].constr.dB)
        {
            return true;
        }
    }

    return false;
}

// ILVES
std::pair<real, real> Constraints::Impl::compute_errors(const t_pbc *const pbc,
                                            ArrayRef<RVec> xprime) {

    // Initialize the largest relative error
    real max_rel_error = 0;
    real avg_rel_error = 0;

    const auto &iatoms = idef->il[F_CONSTR].iatoms;

    // Loop over the n constraints
    size_t n = idef->il[F_CONSTR].size() / 3;
    for (size_t row = 0; row < n; ++row) {
        // Isolate the atoms which partake in bond row
        const int type = iatoms[row * 3 + 0];
        const int a = iatoms[row * 3 + 1];
        const int b = iatoms[row * 3 + 2];

        // Isolate the square of the length of the ith bond
        real sigma2 = gmx::square(idef->iparams[type].constr.dA);

        // Compute the vector from atom a to atom b
        rvec rab;
        if (pbc) {
            pbc_dx(pbc, xprime[b], xprime[a], rab);
        }
        else {
            rvec_sub(xprime[b], xprime[a], rab);
        }

        // Compute the square of the length of r
        const auto scalar = iprod(rab, rab);

        // Compute the relative error (approximation using squares).
        real rel_error = std::abs(0.5 * (sigma2 - scalar) / sigma2);

        // Update the max constraint violation
        max_rel_error = std::max(max_rel_error, rel_error);
        avg_rel_error += rel_error;
    }

    avg_rel_error /= n;

    return std::make_pair(max_rel_error, avg_rel_error);
}

void Constraints::Impl::print_constr_detailed_durations() {
    if (constr_detailed_durations.size() == 0) {
        return;
    }

    std::cout << "Constraint algorithm execution time (ms) per function and thread.\n";

    for (int thread = 0; thread < constr_detailed_durations.size(); ++thread) {
        std::cout << "Thread " << thread << " ==>\n";
        for (auto &duration : constr_detailed_durations[thread]) {
            std::cout << "\t" << duration.first << ": "
                    << duration.second.count() << " ms\n";
        }
    }

    std::cout << "Normalized execution time per function and thread.\n";

    constexpr int setw_size_index = 20;
    constexpr int setw_size_value = 10;

    std::cout << std::setw(setw_size_index) << "";
    for (int thread = 0; thread < constr_detailed_durations.size(); ++thread) {
        std::cout << std::setw(setw_size_value)
                << "Thread " + std::to_string(thread);
    }
    std::cout << "\n";

    for (const auto &iter : constr_detailed_durations[0]) {
        const auto &k = iter.first;
        std::cout << std::setw(setw_size_index) << k;

        double max = 0;
        for (int thread = 0; thread < constr_detailed_durations.size(); ++thread) {
            max = std::max(max, constr_detailed_durations[thread][k].count());
        }
        for (int thread = 0; thread < constr_detailed_durations.size(); ++thread) {
            double val = constr_detailed_durations[thread][k].count() / max;
            std::cout << std::setw(setw_size_value) << std::setprecision(2) << val;
        }
        std::cout << "\n";
    }
}

void Constraints::Impl::mpi_start_printing() {
#ifdef GMX_MPI
    MPI_Comm comm = cr->dd->mpi_comm_all;

    int nranks;
    MPI_Comm_size(comm, &nranks);

    if (nranks < 2) {
        return;
    }

    int rank;
    MPI_Comm_rank(comm, &rank);

    for (int i = 0; i < rank; ++i) {
        int temp = i;
        MPI_Bcast(&temp, 1, MPI_INT, i, comm);
    }

    printf("--- RANK %d ---\n", rank);
#endif
}

void Constraints::Impl::mpi_finish_printing() {
#ifdef GMX_MPI
    MPI_Comm comm = cr->dd->mpi_comm_all;

    int nranks;
    MPI_Comm_size(comm, &nranks);

    if (nranks < 2) {
        return;
    }

    int rank;
    MPI_Comm_rank(comm, &rank);

    for (int i = rank; i < nranks; ++i) {
        int temp = i;
        MPI_Bcast(&temp, 1, MPI_INT, i, comm);
    }
#endif
}
// END ILVES
//! Clears constraint quantities for atoms in nonlocal region.
static void clear_constraint_quantity_nonlocal(gmx_domdec_t* dd, ArrayRef<RVec> q)
{
    int nonlocal_at_start, nonlocal_at_end, at;

    dd_get_constraint_range(dd, &nonlocal_at_start, &nonlocal_at_end);

    for (at = nonlocal_at_start; at < nonlocal_at_end; at++)
    {
        clear_rvec(q[at]);
    }
}

void too_many_constraint_warnings(int eConstrAlg, int warncount)
{
    gmx_fatal(
            FARGS,
            "Too many %s warnings (%d)\n"
            "If you know what you are doing you can %s"
            "set the environment variable GMX_MAXCONSTRWARN to -1,\n"
            "but normally it is better to fix the problem",
            (eConstrAlg == econtLINCS) ? "LINCS" : "SETTLE", warncount,
            (eConstrAlg == econtLINCS) ? "adjust the lincs warning threshold in your mdp file\nor " : "\n");
}

//! Writes out coordinates.
static void write_constr_pdb(const char*          fn,
                             const char*          title,
                             const gmx_mtop_t&    mtop,
                             int                  start,
                             int                  homenr,
                             const t_commrec*     cr,
                             ArrayRef<const RVec> x,
                             const matrix         box)
{
    char          fname[STRLEN];
    FILE*         out;
    int           dd_ac0 = 0, dd_ac1 = 0, i, ii, resnr;
    gmx_domdec_t* dd;
    const char *  anm, *resnm;

    dd = nullptr;
    if (DOMAINDECOMP(cr))
    {
        dd = cr->dd;
        dd_get_constraint_range(dd, &dd_ac0, &dd_ac1);
        start  = 0;
        homenr = dd_ac1;
    }

    if (PAR(cr))
    {
        sprintf(fname, "%s_n%d.pdb", fn, cr->sim_nodeid);
    }
    else
    {
        sprintf(fname, "%s.pdb", fn);
    }

    out = gmx_fio_fopen(fname, "w");

    fprintf(out, "TITLE     %s\n", title);
    gmx_write_pdb_box(out, PbcType::Unset, box);
    int molb = 0;
    for (i = start; i < start + homenr; i++)
    {
        if (dd != nullptr)
        {
            if (i >= dd_numHomeAtoms(*dd) && i < dd_ac0)
            {
                continue;
            }
            ii = dd->globalAtomIndices[i];
        }
        else
        {
            ii = i;
        }
        mtopGetAtomAndResidueName(mtop, ii, &molb, &anm, &resnr, &resnm, nullptr);
        gmx_fprintf_pdb_atomline(out, epdbATOM, ii + 1, anm, ' ', resnm, ' ', resnr, ' ',
                                 10 * x[i][XX], 10 * x[i][YY], 10 * x[i][ZZ], 1.0, 0.0, "");
    }
    fprintf(out, "TER\n");

    gmx_fio_fclose(out);
}

//! Writes out domain contents to help diagnose crashes.
static void dump_confs(FILE*                log,
                       int64_t              step,
                       const gmx_mtop_t&    mtop,
                       int                  start,
                       int                  homenr,
                       const t_commrec*     cr,
                       ArrayRef<const RVec> x,
                       ArrayRef<const RVec> xprime,
                       const matrix         box)
{
    char buf[STRLEN], buf2[22];

    char* env = getenv("GMX_SUPPRESS_DUMP");
    if (env)
    {
        return;
    }

    sprintf(buf, "step%sb", gmx_step_str(step, buf2));
    write_constr_pdb(buf, "initial coordinates", mtop, start, homenr, cr, x, box);
    sprintf(buf, "step%sc", gmx_step_str(step, buf2));
    write_constr_pdb(buf, "coordinates after constraining", mtop, start, homenr, cr, xprime, box);
    if (log)
    {
        fprintf(log, "Wrote pdb files with previous and current coordinates\n");
    }
    fprintf(stderr, "Wrote pdb files with previous and current coordinates\n");
}

bool Constraints::apply(bool                      bLog,
                        bool                      bEner,
                        int64_t                   step,
                        int                       delta_step,
                        real                      step_scaling,
                        ArrayRefWithPadding<RVec> x,
                        ArrayRefWithPadding<RVec> xprime,
                        ArrayRef<RVec>            min_proj,
                        const matrix              box,
                        real                      lambda,
                        real*                     dvdlambda,
                        ArrayRefWithPadding<RVec> v,
                        bool                      computeVirial,
                        tensor                    constraintsVirial,
                        ConstraintVariable        econq)
{
    return impl_->apply(bLog, bEner, step, delta_step, step_scaling, std::move(x),
                        std::move(xprime), min_proj, box, lambda, dvdlambda, std::move(v),
                        computeVirial, constraintsVirial, econq);
}

bool Constraints::Impl::apply(bool                      bLog,
                              bool                      bEner,
                              int64_t                   step,
                              int                       delta_step,
                              real                      step_scaling,
                              ArrayRefWithPadding<RVec> x,
                              ArrayRefWithPadding<RVec> xprime,
                              ArrayRef<RVec>            min_proj,
                              const matrix              box,
                              real                      lambda,
                              real*                     dvdlambda,
                              ArrayRefWithPadding<RVec> v,
                              bool                      computeVirial,
                              tensor                    constraintsVirial,
                              ConstraintVariable        econq)
{
    bool  bOK, bDump;
    int   start;
    real  scaled_delta_t;
    real  invdt, vir_fac = 0, t;
    int   nsettle;
    t_pbc pbc, *pbc_null;
    char  buf[22];
    int   nth;

    wallcycle_start(wcycle, ewcCONSTR);

    if (econq == ConstraintVariable::ForceDispl && !EI_ENERGY_MINIMIZATION(ir.eI))
    {
        gmx_incons(
                "constrain called for forces displacements while not doing energy minimization, "
                "can not do this while the LINCS and SETTLE constraint connection matrices are "
                "mass weighted");
    }

    bOK   = TRUE;
    bDump = FALSE;

    start = 0;

    scaled_delta_t = step_scaling * ir.delta_t;

    /* Prepare time step for use in constraint implementations, and
       avoid generating inf when ir.delta_t = 0. */
    if (ir.delta_t == 0)
    {
        invdt = 0.0;
    }
    else
    {
        invdt = 1.0 / scaled_delta_t;
    }

    if (ir.efep != efepNO && EI_DYNAMICS(ir.eI))
    {
        /* Set the constraint lengths for the step at which this configuration
         * is meant to be. The invmasses should not be changed.
         */
        lambda += delta_step * ir.fepvals->delta_lambda;
    }

    if (computeVirial)
    {
        clear_mat(constraintsVirial);
    }
    const InteractionList& settle = idef->il[F_SETTLE];
    nsettle                       = settle.size() / (1 + NRAL(F_SETTLE));

    if (nsettle > 0)
    {
        nth = gmx_omp_nthreads_get(emntSETTLE);
    }
    else
    {
        nth = 1;
    }

    /* We do not need full pbc when constraints do not cross update groups
     * i.e. when dd->constraint_comm==NULL.
     * Note that PBC for constraints is different from PBC for bondeds.
     * For constraints there is both forward and backward communication.
     */
    if (ir.pbcType != PbcType::No && (cr->dd || pbcHandlingRequired_)
        && !(cr->dd && cr->dd->constraint_comm == nullptr))
    {
        /* With pbc=screw the screw has been changed to a shift
         * by the constraint coordinate communication routine,
         * so that here we can use normal pbc.
         */
        pbc_null = set_pbc_dd(&pbc, ir.pbcType, DOMAINDECOMP(cr) ? cr->dd->numCells : nullptr, FALSE, box);
    }
    else
    {
        pbc_null = nullptr;
    }

    /* Communicate the coordinates required for the non-local constraints
     * for LINCS and/or SETTLE.
     */
    if (cr->dd)
    {
        dd_move_x_constraints(cr->dd, box, x.unpaddedArrayRef(), xprime.unpaddedArrayRef(),
                              econq == ConstraintVariable::Positions);

        if (!v.empty())
        {
            /* We need to initialize the non-local components of v.
             * We never actually use these values, but we do increment them,
             * so we should avoid uninitialized variables and overflows.
             */
            clear_constraint_quantity_nonlocal(cr->dd, v.unpaddedArrayRef());
        }
    }

    // ILVES
    if (CONSTRAINTS_VERBOSE && step % NSCONSTRAINTS_OUTPUT == 0) {
        if (MASTER(cr)) {
            printf("Imposing constraints in time-step: %ld\n", step);
        }
    }
    // END ILVES

    // ILVES
    if (lincsd != nullptr)
    {
        constr_total_calls += 1;

        if (CONSTRAINTS_VERBOSE && step % NSCONSTRAINTS_OUTPUT == 0 && MASTER(cr)) {
            printf("Executing LINCS\n");
        }

        auto start = std::chrono::system_clock::now();
        int niters = 0;
        std::tie(bOK, niters) = constrain_lincs(bLog || bEner,
                                                ir,
                                                step,
                                                lincsd,
                                                inverseMasses_,
                                                ir.lincs_tol,
                                                cr,
                                                ms,
                                                x,
                                                xprime,
                                                min_proj,
                                                box,
                                                pbc_null,
                                                hasMassPerturbedAtoms_,
                                                lambda,
                                                dvdlambda,
                                                invdt,
                                                v.unpaddedArrayRef(),
                                                computeVirial,
                                                constraintsVirial,
                                                econq,
                                                nrnb,
                                                maxwarn,
                                                &warncount_lincs);
        auto end = std::chrono::system_clock::now();
        constr_duration += end - start;

        constr_total_iters += niters;
        if (!bOK && maxwarn < INT_MAX)
        {
            if (log != nullptr)
            {
                fprintf(log, "Constraint error in algorithm %s at step %s\n",
                        econstr_names[econtLINCS], gmx_step_str(step, buf));
            }
            bDump = TRUE;
        }
    }
    // END ILVES

    // ILVES
    if (shaked != nullptr)
    {
        constr_total_calls += 1;

        if (CONSTRAINTS_VERBOSE && step % NSCONSTRAINTS_OUTPUT == 0 && MASTER(cr)) {
            printf("Executing SHAKE\n");
        }

        auto start = std::chrono::system_clock::now();
        int niters = 0;
        std::tie(bOK, niters) = constrain_shake(log,
                                                shaked.get(),
                                                inverseMasses_,
                                                *idef,
                                                ir,
                                                x.unpaddedArrayRef(),
                                                xprime.unpaddedArrayRef(),
                                                min_proj,
                                                pbc_null,
                                                nrnb,
                                                lambda,
                                                dvdlambda,
                                                invdt,
                                                v.unpaddedArrayRef(),
                                                computeVirial,
                                                constraintsVirial,
                                                maxwarn < INT_MAX,
                                                econq);
        auto end = std::chrono::system_clock::now();
        constr_duration += end - start;

        constr_total_iters += niters;

        if (!bOK && maxwarn < INT_MAX) {
            if (log != nullptr) {
                fprintf(log,
                        "Constraint error in algorithm %s at step %s\n",
                        econstr_names[econtSHAKE],
                        gmx_step_str(step, buf));
            }
            bDump = TRUE;
        }
    }
    // END ILVES

    // ILVES
    if (ilves != nullptr) {
        constr_total_calls += 1;

        if (econq != ConstraintVariable::Positions) {
            gmx_fatal(FARGS,
                        "Internal error, ILVES called for constraining something else than "
                        "positions");
        }

        if (ir.efep != efepNO) {
            gmx_fatal(FARGS,
                        "Internal error, ILVES does not support FEP");
        }

        if (CONSTRAINTS_VERBOSE && step % NSCONSTRAINTS_OUTPUT == 0 && MASTER(cr)) {
            printf("Executing ILVES or ILVES-FAST - SIMD %s\n",
                    GMX_SIMD_HAVE_REAL ? "ENABLED" : "DISABLED");
        }

        auto start = std::chrono::system_clock::now();

        int niters = 0;
        std::tie(bOK, niters) = ilves->solve(x.unpaddedArrayRef(),
                                             xprime.unpaddedArrayRef(),
                                             v.unpaddedArrayRef(),
                                             ir.ilves_tol,
                                             scaled_delta_t,
                                             computeVirial,
                                             constraintsVirial,
                                             !v.empty(),
                                             pbc_null,
                                             ir.efep != efepNO,
                                             lambda,
                                             *dvdlambda,
                                             box,
                                             nrnb);

        auto end = std::chrono::system_clock::now();
        constr_duration += end - start;

        constr_total_iters += niters;

        if (!bOK && maxwarn < INT_MAX) {
            if (log != nullptr) {
                fprintf(log,
                        "Constraint error in algorithm ILVES at step %s\n",
                        gmx_step_str(step, buf));
            }
            bDump = TRUE;
        }
    }
    // END ILVES

    // ILVES
    if (CONSTRAINTS_VERBOSE && step % NSCONSTRAINTS_OUTPUT == 0) {
        std::pair<real, real> errors = compute_errors(pbc_null,
                                                        xprime.unpaddedArrayRef());

        if (cr->dd) {
            mpi_start_printing();
        }

        printf("    MAX relative error (approximation) = %E\n", errors.first);
        printf("    AVG relative error (approximation) = %E\n", errors.second);

        if (computeVirial) {
            printf("    Virial: [{%4.6f, %4.6f, %4.6f}, "
                    "{%4.6f, %4.6f, %4.6f}, "
                    "{%4.6f, %4.6f, %4.6f}]\n",
                    constraintsVirial[0][0], constraintsVirial[0][1], constraintsVirial[0][2],
                    constraintsVirial[1][0], constraintsVirial[1][1], constraintsVirial[1][2],
                    constraintsVirial[2][0], constraintsVirial[2][1], constraintsVirial[2][2]);
        }

        // Free energy
        if (ir.efep != efepNO) {
            printf("    Free Energy = %3.15lf\n", *dvdlambda);
        }

        if (cr->dd) {
            mpi_finish_printing();
        }
    }
    // END ILVES

    if (nsettle > 0)
    {
        bool bSettleErrorHasOccurred0 = false;

        switch (econq)
        {
            case ConstraintVariable::Positions:
#pragma omp parallel for num_threads(nth) schedule(static)
                for (int th = 0; th < nth; th++)
                {
                    try
                    {
                        if (th > 0)
                        {
                            clear_mat(threadConstraintsVirial[th]);
                        }

                        csettle(*settled, nth, th, pbc_null, x, xprime, invdt, v, computeVirial,
                                th == 0 ? constraintsVirial : threadConstraintsVirial[th],
                                th == 0 ? &bSettleErrorHasOccurred0 : &bSettleErrorHasOccurred[th]);
                    }
                    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
                }
                inc_nrnb(nrnb, eNR_SETTLE, nsettle);
                if (!v.empty())
                {
                    inc_nrnb(nrnb, eNR_CONSTR_V, nsettle * 3);
                }
                if (computeVirial)
                {
                    inc_nrnb(nrnb, eNR_CONSTR_VIR, nsettle * 3);
                }
                break;
            case ConstraintVariable::Velocities:
            case ConstraintVariable::Derivative:
            case ConstraintVariable::Force:
            case ConstraintVariable::ForceDispl:
#pragma omp parallel for num_threads(nth) schedule(static)
                for (int th = 0; th < nth; th++)
                {
                    try
                    {
                        int calcvir_atom_end;

                        if (!computeVirial)
                        {
                            calcvir_atom_end = 0;
                        }
                        else
                        {
                            calcvir_atom_end = numHomeAtoms_;
                        }

                        if (th > 0)
                        {
                            clear_mat(threadConstraintsVirial[th]);
                        }

                        int start_th = (nsettle * th) / nth;
                        int end_th   = (nsettle * (th + 1)) / nth;

                        if (start_th >= 0 && end_th - start_th > 0)
                        {
                            settle_proj(*settled, econq, end_th - start_th,
                                        settle.iatoms.data() + start_th * (1 + NRAL(F_SETTLE)),
                                        pbc_null, x.unpaddedArrayRef(), xprime.unpaddedArrayRef(),
                                        min_proj, calcvir_atom_end,
                                        th == 0 ? constraintsVirial : threadConstraintsVirial[th]);
                        }
                    }
                    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
                }
                /* This is an overestimate */
                inc_nrnb(nrnb, eNR_SETTLE, nsettle);
                break;
            case ConstraintVariable::Deriv_FlexCon:
                /* Nothing to do, since the are no flexible constraints in settles */
                break;
            default: gmx_incons("Unknown constraint quantity for settle");
        }

        if (computeVirial)
        {
            /* Reduce the virial contributions over the threads */
            for (int th = 1; th < nth; th++)
            {
                m_add(constraintsVirial, threadConstraintsVirial[th], constraintsVirial);
            }
        }

        if (econq == ConstraintVariable::Positions)
        {
            for (int th = 1; th < nth; th++)
            {
                bSettleErrorHasOccurred0 = bSettleErrorHasOccurred0 || bSettleErrorHasOccurred[th];
            }

            if (bSettleErrorHasOccurred0)
            {
                char buf[STRLEN];
                sprintf(buf,
                        "\nstep "
                        "%" PRId64
                        ": One or more water molecules can not be settled.\n"
                        "Check for bad contacts and/or reduce the timestep if appropriate.\n",
                        step);
                if (log)
                {
                    fprintf(log, "%s", buf);
                }
                fprintf(stderr, "%s", buf);
                warncount_settle++;
                if (warncount_settle > maxwarn)
                {
                    too_many_constraint_warnings(-1, warncount_settle);
                }
                bDump = TRUE;

                bOK = FALSE;
            }
        }
    }

    if (computeVirial)
    {
        /* The normal uses of constrain() pass step_scaling = 1.0.
         * The call to constrain() for SD1 that passes step_scaling =
         * 0.5 also passes vir = NULL, so cannot reach this
         * assertion. This assertion should remain until someone knows
         * that this path works for their intended purpose, and then
         * they can use scaled_delta_t instead of ir.delta_t
         * below. */
        assert(gmx_within_tol(step_scaling, 1.0, GMX_REAL_EPS));
        switch (econq)
        {
            case ConstraintVariable::Positions: vir_fac = 0.5 / (ir.delta_t * ir.delta_t); break;
            case ConstraintVariable::Velocities: vir_fac = 0.5 / ir.delta_t; break;
            case ConstraintVariable::Force:
            case ConstraintVariable::ForceDispl: vir_fac = 0.5; break;
            default: gmx_incons("Unsupported constraint quantity for virial");
        }

        if (EI_VV(ir.eI))
        {
            vir_fac *= 2; /* only constraining over half the distance here */
        }
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                constraintsVirial[i][j] *= vir_fac;
            }
        }
    }

    if (bDump)
    {
        dump_confs(log, step, mtop, start, numHomeAtoms_, cr, x.unpaddedArrayRef(),
                   xprime.unpaddedArrayRef(), box);
    }

    if (econq == ConstraintVariable::Positions)
    {
        if (ir.bPull && pull_have_constraint(*pull_work))
        {
            if (EI_DYNAMICS(ir.eI))
            {
                t = ir.init_t + (step + delta_step) * ir.delta_t;
            }
            else
            {
                t = ir.init_t;
            }
            set_pbc(&pbc, ir.pbcType, box);
            pull_constraint(pull_work, masses_, &pbc, cr, ir.delta_t, t,
                            as_rvec_array(x.unpaddedArrayRef().data()),
                            as_rvec_array(xprime.unpaddedArrayRef().data()),
                            as_rvec_array(v.unpaddedArrayRef().data()), constraintsVirial);
        }
        if (ed && delta_step > 0)
        {
            /* apply the essential dynamics constraints here */
            do_edsam(&ir, step, cr, as_rvec_array(xprime.unpaddedArrayRef().data()),
                     as_rvec_array(v.unpaddedArrayRef().data()), box, ed);
        }
    }
    wallcycle_stop(wcycle, ewcCONSTR);

    if (!v.empty() && cFREEZE_)
    {
        /* Set the velocities of frozen dimensions to zero */
        ArrayRef<RVec> vRef = v.unpaddedArrayRef();

        int gmx_unused numThreads = gmx_omp_nthreads_get(emntUpdate);

#pragma omp parallel for num_threads(numThreads) schedule(static)
        for (int i = 0; i < numHomeAtoms_; i++)
        {
            int freezeGroup = cFREEZE_[i];

            for (int d = 0; d < DIM; d++)
            {
                if (ir.opts.nFreeze[freezeGroup][d])
                {
                    vRef[i][d] = 0;
                }
            }
        }
    }

    return bOK;
} // namespace gmx

ArrayRef<real> Constraints::rmsdData() const
{
    if (impl_->lincsd)
    {
        return lincs_rmsdData(impl_->lincsd);
    }
    else
    {
        return {};
    }
}

real Constraints::rmsd() const
{
    if (impl_->lincsd)
    {
        return lincs_rmsd(impl_->lincsd);
    }
    else
    {
        return 0;
    }
}

int Constraints::numConstraintsTotal()
{
    return impl_->ncon_tot;
}

FlexibleConstraintTreatment flexibleConstraintTreatment(bool haveDynamicsIntegrator)
{
    if (haveDynamicsIntegrator)
    {
        return FlexibleConstraintTreatment::Include;
    }
    else
    {
        return FlexibleConstraintTreatment::Exclude;
    }
}

/*! \brief Returns a block struct to go from atoms to constraints
 *
 * The block struct will contain constraint indices with lower indices
 * directly matching the order in F_CONSTR and higher indices matching
 * the order in F_CONSTRNC offset by the number of constraints in F_CONSTR.
 *
 * \param[in]  numAtoms  The number of atoms to construct the list for
 * \param[in]  ilists    The interaction lists, size F_NRE
 * \param[in]  iparams   Interaction parameters, can be null when
 *                       \p flexibleConstraintTreatment==Include
 * \param[in]  flexibleConstraintTreatment  The flexible constraint treatment,
 *                                          see enum above
 *
 * \returns a block struct with all constraints for each atom
 */
static ListOfLists<int> makeAtomsToConstraintsList(int                             numAtoms,
                                                   ArrayRef<const InteractionList> ilists,
                                                   ArrayRef<const t_iparams>       iparams,
                                                   FlexibleConstraintTreatment flexibleConstraintTreatment)
{
    GMX_ASSERT(flexibleConstraintTreatment == FlexibleConstraintTreatment::Include || !iparams.empty(),
               "With flexible constraint detection we need valid iparams");

    std::vector<int> count(numAtoms);

    for (int ftype = F_CONSTR; ftype <= F_CONSTRNC; ftype++)
    {
        const InteractionList& ilist  = ilists[ftype];
        const int              stride = 1 + NRAL(ftype);
        for (int i = 0; i < ilist.size(); i += stride)
        {
            if (flexibleConstraintTreatment == FlexibleConstraintTreatment::Include
                || !isConstraintFlexible(iparams, ilist.iatoms[i]))
            {
                for (int j = 1; j < 3; j++)
                {
                    int a = ilist.iatoms[i + j];
                    count[a]++;
                }
            }
        }
    }

    std::vector<int> listRanges(numAtoms + 1);
    for (int a = 0; a < numAtoms; a++)
    {
        listRanges[a + 1] = listRanges[a] + count[a];
        count[a]          = 0;
    }
    std::vector<int> elements(listRanges[numAtoms]);

    /* The F_CONSTRNC constraints have constraint numbers
     * that continue after the last F_CONSTR constraint.
     */
    int numConstraints = 0;
    for (int ftype = F_CONSTR; ftype <= F_CONSTRNC; ftype++)
    {
        const InteractionList& ilist  = ilists[ftype];
        const int              stride = 1 + NRAL(ftype);
        for (int i = 0; i < ilist.size(); i += stride)
        {
            if (flexibleConstraintTreatment == FlexibleConstraintTreatment::Include
                || !isConstraintFlexible(iparams, ilist.iatoms[i]))
            {
                for (int j = 1; j < 3; j++)
                {
                    const int a                          = ilist.iatoms[i + j];
                    elements[listRanges[a] + count[a]++] = numConstraints;
                }
            }
            numConstraints++;
        }
    }

    return ListOfLists<int>(std::move(listRanges), std::move(elements));
}

ListOfLists<int> make_at2con(int                             numAtoms,
                             ArrayRef<const InteractionList> ilist,
                             ArrayRef<const t_iparams>       iparams,
                             FlexibleConstraintTreatment     flexibleConstraintTreatment)
{
    return makeAtomsToConstraintsList(numAtoms, ilist, iparams, flexibleConstraintTreatment);
}

ListOfLists<int> make_at2con(const gmx_moltype_t&           moltype,
                             gmx::ArrayRef<const t_iparams> iparams,
                             FlexibleConstraintTreatment    flexibleConstraintTreatment)
{
    return makeAtomsToConstraintsList(moltype.atoms.nr, makeConstArrayRef(moltype.ilist), iparams,
                                      flexibleConstraintTreatment);
}

//! Return the number of flexible constraints in the \c ilist and \c iparams.
int countFlexibleConstraints(ArrayRef<const InteractionList> ilist, ArrayRef<const t_iparams> iparams)
{
    int nflexcon = 0;
    for (int ftype = F_CONSTR; ftype <= F_CONSTRNC; ftype++)
    {
        const int numIatomsPerConstraint = 3;
        for (int i = 0; i < ilist[ftype].size(); i += numIatomsPerConstraint)
        {
            const int type = ilist[ftype].iatoms[i];
            if (iparams[type].constr.dA == 0 && iparams[type].constr.dB == 0)
            {
                nflexcon++;
            }
        }
    }

    return nflexcon;
}

//! Returns the index of the settle to which each atom belongs.
static std::vector<int> make_at2settle(int natoms, const InteractionList& ilist)
{
    /* Set all to no settle */
    std::vector<int> at2s(natoms, -1);

    const int stride = 1 + NRAL(F_SETTLE);

    for (int s = 0; s < ilist.size(); s += stride)
    {
        at2s[ilist.iatoms[s + 1]] = s / stride;
        at2s[ilist.iatoms[s + 2]] = s / stride;
        at2s[ilist.iatoms[s + 3]] = s / stride;
    }

    return at2s;
}

void Constraints::Impl::setConstraints(gmx_localtop_t* top,
                                       int             numAtoms,
                                       int             numHomeAtoms,
                                       real*           masses,
                                       real*           inverseMasses,
                                       const bool      hasMassPerturbedAtoms,
                                       const real      lambda,
                                       unsigned short* cFREEZE)
{
    numAtoms_              = numAtoms;
    numHomeAtoms_          = numHomeAtoms;
    masses_                = masses;
    inverseMasses_         = inverseMasses;
    hasMassPerturbedAtoms_ = hasMassPerturbedAtoms;
    lambda_                = lambda;
    cFREEZE_               = cFREEZE;

    idef = &top->idef;

    if (ncon_tot > 0)
    {
        // ILVES
        auto start = std::chrono::system_clock::now();

        if (ir.eConstrAlg == econtILVES) {
            ilves = std::make_unique<IlvesAsym>(cr,
                                                *idef,
                                                inverseMasses_,
                                                omp_get_max_threads(),
                                                constr_detailed_durations);
        }
        if (ir.eConstrAlg == econtILVESF) {
            ilves = std::make_unique<IlvesSym>(cr, 
                                                *idef,
                                                inverseMasses_,
                                                omp_get_max_threads(),
                                                constr_detailed_durations);
        }
        /* With DD we might also need to call LINCS on a domain no constraints for
         * communicating coordinates to other nodes that do have constraints.
         */
        if (ir.eConstrAlg == econtLINCS)
        {
            set_lincs(*idef, numAtoms_, inverseMasses_, lambda_, EI_DYNAMICS(ir.eI), cr, lincsd);
        }
        if (ir.eConstrAlg == econtSHAKE)
        {
            if (cr->dd)
            {
                // We are using the local topology, so there are only
                // F_CONSTR constraints.
                GMX_RELEASE_ASSERT(idef->il[F_CONSTRNC].empty(),
                                   "Here we should not have no-connect constraints");
                make_shake_sblock_dd(shaked.get(), idef->il[F_CONSTR], CONSTRAINTS_VERBOSE == 2);
            }
            else
            {
                make_shake_sblock_serial(shaked.get(), &top->idef, numAtoms_, CONSTRAINTS_VERBOSE == 2);
            }
        }

        auto end = std::chrono::system_clock::now();
        constr_init_duration += end - start;
        // END ILVES
    }

    if (settled)
    {
        settled->setConstraints(idef->il[F_SETTLE], numHomeAtoms_, masses_, inverseMasses_);
    }

    /* Make a selection of the local atoms for essential dynamics */
    if (ed && cr->dd)
    {
        dd_make_local_ed_indices(cr->dd, ed);
    }
}

void Constraints::setConstraints(gmx_localtop_t* top,
                                 const int       numAtoms,
                                 const int       numHomeAtoms,
                                 real*           masses,
                                 real*           inverseMasses,
                                 const bool      hasMassPerturbedAtoms,
                                 const real      lambda,
                                 unsigned short* cFREEZE)
{
    impl_->setConstraints(top, numAtoms, numHomeAtoms, masses, inverseMasses, hasMassPerturbedAtoms,
                          lambda, cFREEZE);
}

/*! \brief Makes a per-moleculetype container of mappings from atom
 * indices to constraint indices.
 *
 * Note that flexible constraints are only enabled with a dynamical integrator. */
static std::vector<ListOfLists<int>> makeAtomToConstraintMappings(const gmx_mtop_t& mtop,
                                                                  FlexibleConstraintTreatment flexibleConstraintTreatment)
{
    std::vector<ListOfLists<int>> mapping;
    mapping.reserve(mtop.moltype.size());
    for (const gmx_moltype_t& moltype : mtop.moltype)
    {
        mapping.push_back(make_at2con(moltype, mtop.ffparams.iparams, flexibleConstraintTreatment));
    }
    return mapping;
}

Constraints::Constraints(const gmx_mtop_t&     mtop,
                         const t_inputrec&     ir,
                         pull_t*               pull_work,
                         FILE*                 log,
                         const t_commrec*      cr,
                         const gmx_multisim_t* ms,
                         t_nrnb*               nrnb,
                         gmx_wallcycle*        wcycle,
                         bool                  pbcHandlingRequired,
                         int                   numConstraints,
                         int                   numSettles) :
    impl_(new Impl(mtop, ir, pull_work, log, cr, ms, nrnb, wcycle, pbcHandlingRequired, numConstraints, numSettles))
{
}

Constraints::Impl::Impl(const gmx_mtop_t&     mtop_p,
                        const t_inputrec&     ir_p,
                        pull_t*               pull_work,
                        FILE*                 log_p,
                        const t_commrec*      cr_p,
                        const gmx_multisim_t* ms_p,
                        t_nrnb*               nrnb_p,
                        gmx_wallcycle*        wcycle_p,
                        bool                  pbcHandlingRequired,
                        int                   numConstraints,
                        int                   numSettles) :
    ncon_tot(numConstraints),
    mtop(mtop_p),
    pbcHandlingRequired_(pbcHandlingRequired),
    log(log_p),
    cr(cr_p),
    ms(ms_p),
    pull_work(pull_work),
    ir(ir_p),
    nrnb(nrnb_p),
    wcycle(wcycle_p)
{
    if (numConstraints + numSettles > 0 && ir.epc == epcMTTK)
    {
        gmx_fatal(FARGS, "Constraints are not implemented with MTTK pressure control.");
    }

    // ILVES
    char *CONSTRAINTS_VERBOSE_ev = std::getenv("CONSTRAINTS_VERBOSE");
    if (CONSTRAINTS_VERBOSE_ev == NULL) {
        CONSTRAINTS_VERBOSE = 0;
    }
    else {
        int CONSTRAINTS_VERBOSE_int = atoi(CONSTRAINTS_VERBOSE_ev);
        CONSTRAINTS_VERBOSE = (CONSTRAINTS_VERBOSE_int >= 0 && CONSTRAINTS_VERBOSE_int <= 2)
                                ? CONSTRAINTS_VERBOSE_int : 0;
    }

    char *NSCONSTRAINTS_OUTPUT_ev = std::getenv("NSCONSTRAINTS_OUTPUT");
    if (NSCONSTRAINTS_OUTPUT_ev == NULL) {
        NSCONSTRAINTS_OUTPUT = 1;
    }
    else {
        int NSCONSTRAINTS_OUTPUT_int = atoi(NSCONSTRAINTS_OUTPUT_ev);
        NSCONSTRAINTS_OUTPUT = (NSCONSTRAINTS_OUTPUT_int >= 1) ?
                                NSCONSTRAINTS_OUTPUT_int : 1;
    }

    // END ILVES

    nflexcon = 0;
    if (numConstraints > 0)
    {
        at2con_mt = makeAtomToConstraintMappings(mtop, flexibleConstraintTreatment(EI_DYNAMICS(ir.eI)));

        for (const gmx_molblock_t& molblock : mtop.molblock)
        {
            int count = countFlexibleConstraints(mtop.moltype[molblock.type].ilist, mtop.ffparams.iparams);
            nflexcon += molblock.nmol * count;
        }

        if (nflexcon > 0)
        {
            if (log)
            {
                fprintf(log, "There are %d flexible constraints\n", nflexcon);
                if (ir.fc_stepsize == 0)
                {
                    fprintf(log,
                            "\n"
                            "WARNING: step size for flexible constraining = 0\n"
                            "         All flexible constraints will be rigid.\n"
                            "         Will try to keep all flexible constraints at their original "
                            "length,\n"
                            "         but the lengths may exhibit some drift.\n\n");
                    nflexcon = 0;
                }
            }
            if (nflexcon > 0)
            {
                please_cite(log, "Hess2002");
            }
        }

        if (ir.eConstrAlg == econtLINCS)
        {
            lincsd = init_lincs(log, mtop, nflexcon, at2con_mt,
                                DOMAINDECOMP(cr) && ddHaveSplitConstraints(*cr->dd), ir.nLincsIter,
                                ir.nProjOrder);
        }

        if (ir.eConstrAlg == econtSHAKE)
        {
            if (DOMAINDECOMP(cr) && ddHaveSplitConstraints(*cr->dd))
            {
                gmx_fatal(FARGS,
                          "SHAKE is not supported with domain decomposition and constraint that "
                          "cross domain boundaries, use LINCS");
            }
            if (nflexcon)
            {
                gmx_fatal(FARGS,
                          "For this system also velocities and/or forces need to be constrained, "
                          "this can not be done with SHAKE, you should select LINCS");
            }
            please_cite(log, "Ryckaert77a");
            if (ir.bShakeSOR)
            {
                please_cite(log, "Barth95a");
            }

            shaked = std::make_unique<shakedata>();
        }
    }

    if (numSettles > 0)
    {
        please_cite(log, "Miyamoto92a");

        settled = std::make_unique<SettleData>(mtop);

        /* Make an atom to settle index for use in domain decomposition */
        for (size_t mt = 0; mt < mtop.moltype.size(); mt++)
        {
            at2settle_mt.emplace_back(
                    make_at2settle(mtop.moltype[mt].atoms.nr, mtop.moltype[mt].ilist[F_SETTLE]));
        }

        /* Allocate thread-local work arrays */
        int nthreads = gmx_omp_nthreads_get(emntSETTLE);
        if (nthreads > 1 && threadConstraintsVirial == nullptr)
        {
            snew(threadConstraintsVirial, nthreads);
            snew(bSettleErrorHasOccurred, nthreads);
        }
    }

    maxwarn   = 999;
    char* env = getenv("GMX_MAXCONSTRWARN");
    if (env)
    {
        maxwarn = 0;
        sscanf(env, "%8d", &maxwarn);
        if (maxwarn < 0)
        {
            maxwarn = INT_MAX;
        }
        if (log)
        {
            fprintf(log, "Setting the maximum number of constraint warnings to %d\n", maxwarn);
        }
        if (MASTER(cr))
        {
            fprintf(stderr, "Setting the maximum number of constraint warnings to %d\n", maxwarn);
        }
    }
    warncount_lincs  = 0;
    warncount_settle = 0;
}

Constraints::Impl::~Impl()
{
    if (ncon_tot > 0) {
        // Print in order
        if (cr->dd) {
            mpi_start_printing();
        }
        printf("CONSTR ALGORITHM time (s) = %lf\n", constr_duration.count());
        printf("CONSTR ALGORITHM initialization time (s) = %lf\n", constr_init_duration.count());
        printf("CONSTR ALGORITHM iterations = %d\n", constr_total_iters);
        printf("CONSTR ALGORITHM total calls = %d\n", constr_total_calls);
        printf("CONSTR ALGORITHM average iterations = %lf\n", (double)constr_total_iters / constr_total_calls);
        print_constr_detailed_durations();
        if (cr->dd) {
            mpi_finish_printing();
        }
    }

    if (bSettleErrorHasOccurred != nullptr)
    {
        sfree(bSettleErrorHasOccurred);
    }
    if (threadConstraintsVirial != nullptr)
    {
        sfree(threadConstraintsVirial);
    }
    done_lincs(lincsd);
}

void Constraints::saveEdsamPointer(gmx_edsam* ed)
{
    impl_->ed = ed;
}

ArrayRef<const ListOfLists<int>> Constraints::atom2constraints_moltype() const
{
    return impl_->at2con_mt;
}

ArrayRef<const std::vector<int>> Constraints::atom2settle_moltype() const
{
    return impl_->at2settle_mt;
}

void do_constrain_first(FILE*                     fplog,
                        gmx::Constraints*         constr,
                        const t_inputrec*         ir,
                        int                       numAtoms,
                        int                       numHomeAtoms,
                        ArrayRefWithPadding<RVec> x,
                        ArrayRefWithPadding<RVec> v,
                        const matrix              box,
                        real                      lambda)
{
    int     i, m, start, end;
    int64_t step;
    real    dt = ir->delta_t;
    real    dvdl_dum;

    PaddedVector<RVec> savex(numAtoms);

    start = 0;
    end   = numHomeAtoms;

    if (debug)
    {
        fprintf(debug, "vcm: start=%d, homenr=%d, end=%d\n", start, numHomeAtoms, end);
    }
    /* Do a first constrain to reset particles... */
    step = ir->init_step;
    if (fplog)
    {
        char buf[STEPSTRSIZE];
        fprintf(fplog, "\nConstraining the starting coordinates (step %s)\n", gmx_step_str(step, buf));
    }
    dvdl_dum = 0;

    bool needsLogging  = true;
    bool computeEnergy = false;
    bool computeVirial = false;
    /* constrain the current position */
    constr->apply(needsLogging, computeEnergy, step, 0, 1.0, x, x, {}, box, lambda, &dvdl_dum, {},
                  computeVirial, nullptr, gmx::ConstraintVariable::Positions);
    if (EI_VV(ir->eI))
    {
        /* constrain the inital velocity, and save it */
        /* also may be useful if we need the ekin from the halfstep for velocity verlet */
        constr->apply(needsLogging, computeEnergy, step, 0, 1.0, x, v, v.unpaddedArrayRef(), box, lambda,
                      &dvdl_dum, {}, computeVirial, nullptr, gmx::ConstraintVariable::Velocities);
    }
    /* constrain the inital velocities at t-dt/2 */
    if (EI_STATE_VELOCITY(ir->eI) && ir->eI != eiVV)
    {
        auto subX = x.paddedArrayRef().subArray(start, end);
        auto subV = v.paddedArrayRef().subArray(start, end);
        for (i = start; (i < end); i++)
        {
            for (m = 0; (m < DIM); m++)
            {
                /* Reverse the velocity */
                subV[i][m] = -subV[i][m];
                /* Store the position at t-dt in buf */
                savex[i][m] = subX[i][m] + dt * subV[i][m];
            }
        }
        /* Shake the positions at t=-dt with the positions at t=0
         * as reference coordinates.
         */
        if (fplog)
        {
            char buf[STEPSTRSIZE];
            fprintf(fplog, "\nConstraining the coordinates at t0-dt (step %s)\n", gmx_step_str(step, buf));
        }
        dvdl_dum = 0;
        constr->apply(needsLogging, computeEnergy, step, -1, 1.0, x, savex.arrayRefWithPadding(),
                      {}, box, lambda, &dvdl_dum, v, computeVirial, nullptr,
                      gmx::ConstraintVariable::Positions);

        for (i = start; i < end; i++)
        {
            for (m = 0; m < DIM; m++)
            {
                /* Re-reverse the velocities */
                subV[i][m] = -subV[i][m];
            }
        }
    }
}

void constrain_velocities(gmx::Constraints* constr,
                          bool              do_log,
                          bool              do_ene,
                          int64_t           step,
                          t_state*          state,
                          real*             dvdlambda,
                          gmx_bool          computeVirial,
                          tensor            constraintsVirial)
{
    if (constr != nullptr)
    {
        constr->apply(do_log, do_ene, step, 1, 1.0, state->x.arrayRefWithPadding(),
                      state->v.arrayRefWithPadding(), state->v.arrayRefWithPadding().unpaddedArrayRef(),
                      state->box, state->lambda[efptBONDED], dvdlambda, ArrayRefWithPadding<RVec>(),
                      computeVirial, constraintsVirial, ConstraintVariable::Velocities);
    }
}

void constrain_coordinates(gmx::Constraints*         constr,
                           bool                      do_log,
                           bool                      do_ene,
                           int64_t                   step,
                           t_state*                  state,
                           ArrayRefWithPadding<RVec> xp,
                           real*                     dvdlambda,
                           gmx_bool                  computeVirial,
                           tensor                    constraintsVirial)
{
    if (constr != nullptr)
    {
        constr->apply(do_log, do_ene, step, 1, 1.0, state->x.arrayRefWithPadding(), std::move(xp),
                      ArrayRef<RVec>(), state->box, state->lambda[efptBONDED], dvdlambda,
                      state->v.arrayRefWithPadding(), computeVirial, constraintsVirial,
                      ConstraintVariable::Positions);
    }
}

} // namespace gmx
