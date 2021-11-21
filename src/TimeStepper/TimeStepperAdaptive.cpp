/**
 * Implementation of an adaptive time stepper for DREAM.
 *
 * The stepper advances the system in time, occasionally (or always) checking
 * that the solution is converged w.r.t. the time step dt. This is achieved by
 * feeding the appropriate time step to the solver when 'NextTime()' is called.
 * To determine whether to increase or decrease the time step, the stepper first
 * takes two steps of size 'dt/2', then returns to the initial state and takes
 * a single step of size 'dt'. By comparing the resulting solutions in these
 * two cases, an estimate for the truncation error can be obtained.
 *
 * To keep track of where in the stepping we are, the stepper moves between a
 * number of stages during the stepping. The stages are:
 *
 *   NORMAL       -- Take a normal time step of size 'dt'. No error checks are
 *                   done after such a step.
 *   FIRST_HALF   -- Take the first of the two steps of size 'dt/2'.
 *   SECOND_HALF  -- Take the second of the two steps of size 'dt/2'.
 *   FULL         -- Take a full step of size 'dt' and check the error after
 *                   taking it. Based on the result of the check, we also
 *                   update the time step 'dt'.
 * 
 * The progression of stages can be illustrated as follows:
 *
 *  -> NORMAL -> ... -> FIRST_HALF -> SECOND_HALF -> FULL ---> [converged, increase dt] -----> |
 * |                                                       |                                   |
 * |                                                        -> [not converged, decrease dt] -> |
 * |___________________________________________________________________________________________|
 */

#include <iostream>
#include <vector>
#include "DREAM/IO.hpp"
#include "DREAM/TimeStepper/TimeStepperAdaptive.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 *
 * tMax:         Final simulation time.
 * dt0:          Initial time step.
 * uqh:          UnknownQuantityHandler of solver.
 * nontrivials:  List of non-trivial unknowns.
 * cc:           Object to use for checking time stepper convergence.
 * checkEvery:   Number of time steps to take _without_ doing a convergence
 *               check after each check (i.e. 0 => check _every_ time step,
 *               1 => check every other etc.)
 * constantStep: Overrides the adaptive stepper and ensures that the time
 *               step is held constant. This parameter is solely intended
 *               for debugging and verification of the stepper.
 */
TimeStepperAdaptive::TimeStepperAdaptive(
    const real_t tMax, const real_t dt0, FVM::UnknownQuantityHandler *uqh,
    EquationSystem *eqsys,
    vector<len_t>& nontrivials, ConvergenceChecker *cc, int_t checkEvery,
    bool verbose, bool constantStep
) : TimeStepper(uqh, eqsys), tMax(tMax), dt(dt0), nontrivials(nontrivials),
  checkEvery(checkEvery), verbose(verbose), constantStep(constantStep) {
    
    this->stepsSinceCheck = checkEvery;

    if (cc == nullptr) {
        const real_t RELTOL = 1e-6;
        this->convChecker = new ConvergenceChecker(uqh, nontrivials, RELTOL);
    } else
        this->convChecker = cc;

    // Initial guess is to solve system in a single step...
    if (dt > tMax)
        this->dt = this->tMax;

    if (verbose)
        DREAM::IO::PrintInfo(
            "[TimeStepper] initial dt = %.6e", this->dt
        );
}

/**
 * Destructor.
 */
TimeStepperAdaptive::~TimeStepperAdaptive() {
    DeallocateSolutions();

    delete this->convChecker;
}


/**
 * Advance the solver to the next stage in the time
 * stepping (the sequence is
 *
 *   NORMAL -> ... -> 1st half -> 2nd half -> full -> NORMAL -> ...
 * )
 */
TimeStepperAdaptive::ts_stage TimeStepperAdaptive::AdvanceStage() {
    ts_stage stg = this->currentStage;

    // If an exception occured, last we have rolled back the
    // solution and should retry the step. That means we should
    // behave as if the last step was 'NORMAL' (and we override
    // the 'ShouldCheckError()' in this case)...
    bool overrideCheckError = false;
    if (this->stepsWithException > 0) {
        stg = STAGE_NORMAL;
        overrideCheckError = true;
        this->currentTime = this->initTime;
    }
    
    switch (stg) {
        case STAGE_NORMAL:
            if (ShouldCheckError() || overrideCheckError) {
                this->stepsSinceCheck = 0;
                stg = STAGE_FIRST_HALF;
                this->unknowns->SaveStep(this->initTime, false);
                this->initTime = this->currentTime;
            } else {
                // Keep taking normal steps
                this->stepsSinceCheck++;
                this->currentTime += this->dt;
            }
            break;

        case STAGE_FIRST_HALF:
            stg = STAGE_SECOND_HALF;
            this->currentTime += 0.5*this->dt;
            break;

        case STAGE_SECOND_HALF:
            stg = STAGE_FULL;
            // Get solution due to two half steps
            CopySolution(&this->sol_half);
            RestoreInitialSolution(3);
            this->currentTime = this->initTime;
            break;

        case STAGE_FULL: {
            // This means that the previous step was the full step,
            // so that we should now evaluate how successful the time
            // stepping has been. We must also update the current time
            // (and make sure it's done before modifying the time step dt)
            if (this->stepSucceeded) {
                stg = (checkEvery > 0) ? STAGE_NORMAL : STAGE_FIRST_HALF;
                this->currentTime += this->oldDt;
                this->currentStep++;

                if (stg == STAGE_FIRST_HALF) {
                    this->initTime = this->currentTime;
                    this->unknowns->SaveStep(this->initTime, false);
                }
            } else {  // Time-stepping failed => redo step
                stg = STAGE_FIRST_HALF;
                RestoreInitialSolution(2);
            }
        } break;

        default:
            throw TimeStepperException("Unrecognized stage: %d.");
    }

    if (this->verbose)
        DREAM::IO::PrintInfo(
            "[TimeStepper] Advancing from stage %s -> %s",
            GetStageName(this->currentStage), GetStageName(stg)
        );

    this->currentStage = stg;

    return stg;
}

/**
 * Allocate memory for storing the initial and half-step
 * solution vectors.
 *
 * full_size:    Elements in vector for storing data for _all_ unknown quantities.
 * nontriv_size: Elements in vector for storing data from only non-trivial unknowns.
 */
void TimeStepperAdaptive::AllocateSolutions(const len_t size) {
    if (this->sol_half != nullptr)
        DeallocateSolutions();

    this->sol_size = size;
    this->sol_init = new real_t[size];
    this->sol_half = new real_t[size];
    this->sol_full = new real_t[size];
}

/**
 * Deallocate memory for the solution vectors.
 */
void TimeStepperAdaptive::DeallocateSolutions() {
    if (this->sol_init != nullptr)
        delete [] this->sol_init;
    if (this->sol_half != nullptr)
        delete [] this->sol_half;
    if (this->sol_full != nullptr)
        delete [] this->sol_full;
}

/**
 * Copy solution vector to temporary storage.
 */
void TimeStepperAdaptive::CopySolution(real_t **sol) {
    const len_t SSIZE = this->unknowns->GetLongVectorSize(this->nontrivials);
    if (sol_size != SSIZE) {
        AllocateSolutions(SSIZE);
    }

    this->unknowns->GetLongVector(this->nontrivials, *sol);
}

/**
 * Returns the current time.
 */
real_t TimeStepperAdaptive::CurrentTime() const {
    return this->currentTime;
}

/**
 * Returns the name of the specified time stepper stage.
 *
 * stg: Stage to provide name for.
 */
const char *TimeStepperAdaptive::GetStageName(ts_stage stg) {
    switch (stg) {
        case STAGE_NORMAL: return "NORMAL";
        case STAGE_FIRST_HALF: return "FIRST_HALF";
        case STAGE_SECOND_HALF: return "SECOND_HALF";
        case STAGE_FULL: return "FULL";

        default: return "<INVALID>";
    }
}

/**
 * This method is called when an exception was thrown while
 * the next time step was being taken.
 *
 * In the adaptive time stepper, we rollback to the initial
 * state (unless the exception has been caught too many times)
 * and try again with an even shorter time step.
 *
 * ex: The exception that was caught.
 */
void TimeStepperAdaptive::HandleException(FVM::FVMException &ex) {
    // Maximum number of exceptions thrown?
    if (this->stepsWithException >= MAX_STEPS_WITH_EXCEPTION) {
        DREAM::IO::PrintError("TimeStepper: Caught exception for the last time. Rethrowing...");
        throw ex;
    }

    // Count the exception
    this->stepsWithException++;

    // Determine how many steps to roll back...
    switch (this->currentStage) {
        case STAGE_NORMAL:
            //RestoreInitialSolution(1);
            // Cannot roll back from this state (and it's unlikely that we
            // will ever end up here), so we rethrow instead...
            DREAM::IO::PrintError("TimeStepper: Cannot roll back from stage 'NORMAL'. Rethrowing exception...");
            throw ex;

        case STAGE_FIRST_HALF:
            RestoreInitialSolution(1, false);
            break;

        case STAGE_SECOND_HALF:
            RestoreInitialSolution(2, false);
            break;

        case STAGE_FULL:
            RestoreInitialSolution(1, false);
            break;

        default:
            throw TimeStepperException(
                "TimeStepper::HandleException: Unrecognized stage: " LEN_T_PRINTF_FMT ".",
                this->currentStage
            );
    }

    // Since we have no way of determining how small the time
    // step should be (we just know that the current time step
    // is too large and crashes the solver), we simply multiply  
    // it by some ad-hoc factor <1
    this->dt *= STEP_REDUCTION_AT_EXCEPTION;

    if (this->verbose) {
        DREAM::IO::PrintInfo("Caught exception. Reducing time step to %e", this->dt);
    }
}

/**
 * Returns 'true' if the simulation has reached the final
 * time. Returns 'false' otherwise.
 */
bool TimeStepperAdaptive::IsFinished() {
    bool v = (this->stepSucceeded && (this->currentTime+this->oldDt) >= this->tMax);
#ifdef DREAM_IS_PYTHON_LIBRARY
    return (v || this->PythonIsTerminate());
#else
    return v;
#endif
}

/**
 * Returns 'true' if the result of the current time step
 * should be saved to the output.
 */
bool TimeStepperAdaptive::IsSaveStep() {
    return 
        (this->currentStage == STAGE_NORMAL || this->currentStage == STAGE_FULL) &&
        this->stepSucceeded;
}

/**
 * Returns the maximum simulation time for this simulation.
 */
real_t TimeStepperAdaptive::MaxTime() const {
    return this->tMax;
}

/**
 * Calculate and return the time of the next step to take.
 * The asbolute time value is returned and the caller will
 * have to turn that into a time step themselves.
 */
real_t TimeStepperAdaptive::NextTime() {
    ts_stage stg = AdvanceStage();
    real_t newTime = this->initTime;

    switch (stg) {
        case STAGE_NORMAL:
            newTime += this->dt;
            break;

        case STAGE_FIRST_HALF:
            newTime += 0.5 * this->dt;
            break;

        case STAGE_SECOND_HALF:
            newTime += this->dt;
            break;

        case STAGE_FULL:
            newTime += this->dt;
            break;
    }

    return newTime;
}

/**
 * Print current time stepping progress.
 */
void TimeStepperAdaptive::PrintProgress() {
    if (!this->stepSucceeded)
        return;

    const len_t PERC_FMT_PREC = 2;      // Precision (after decimal point) in percentage
    //                          100 . XX            %
    const len_t PERC_FMT_LENGTH = 3+1+PERC_FMT_PREC+1;
    const len_t EDGE_LENGTH = 1;
    const len_t PROG_LENGTH = PROGRESSBAR_LENGTH-2*EDGE_LENGTH - PERC_FMT_LENGTH - 1;

    cout << "\r[";
    real_t perc     = (CurrentTime() + this->oldDt)/this->tMax;
    len_t threshold = static_cast<len_t>(perc * PROG_LENGTH);
    
    for (len_t i = 0; i < PROG_LENGTH; i++) {
        if (i < threshold)
            cout << '#';
        else
            cout << '-';
    }

    cout << "] ";
    printf(
        "%*.*f%% (step " LEN_T_PRINTF_FMT ", dt = %.5e)",
        int(4+PERC_FMT_PREC), int(PERC_FMT_PREC),
        perc*100.0, this->currentStep, this->dt
    );

    // Ensure that output is written right away (otherwise it may
    // not be written until the end-of-line character is written)
    cout << flush;
}

/**
 * Check whether the last step taken (if it was a full step) reached
 * the desired tolerance.
 */
void TimeStepperAdaptive::ValidateStep() {
    // Reset exception counter (because this method is only called
    // if the solver succeeded)
    this->stepsWithException = 0;

    if (this->currentStage == STAGE_FULL) {
        if (UpdateStep())
            this->stepSucceeded = true;
        else
            this->stepSucceeded = false;
    } else
        this->stepSucceeded = (this->currentStage == STAGE_NORMAL);
}

/**
 * Restore the initial solution.
 *
 * nSteps:   Number of time steps to roll back.
 * pushinit: If 'true', pushes the restored "initial" solution
 *           after rolling back.
 */
void TimeStepperAdaptive::RestoreInitialSolution(const len_t nSteps, bool pushinit) {
    for (len_t i = 0; i < nSteps; i++)
        this->unknowns->RollbackSaveStep();

    if (this->sol_init == nullptr) {
        const len_t SSIZE = this->unknowns->GetLongVectorSize(this->nontrivials);
        AllocateSolutions(SSIZE);
    }

    // Restore initial solution in solver
    this->unknowns->GetLongVector(this->nontrivials, this->sol_init);
    this->solver->SetInitialGuess(this->sol_init);

    // Save current "initial" step
    if (pushinit)
        this->unknowns->SaveStep(this->initTime, false);
}

/**
 * Returns 'true' if the tolerance should be checked before/while
 * taking the next time step.
 */
bool TimeStepperAdaptive::ShouldCheckError() {
    return (
        this->stepsSinceCheck >= this->checkEvery
    );
}

/**
 * Update the size of the time step 'dt' based on
 * the error due to the previous step.
 */
bool TimeStepperAdaptive::UpdateStep() {
    // Store old dt for safekeeping
    this->oldDt = this->dt;

    CopySolution(&this->sol_full);

    bool converged = this->convChecker->IsConverged(this->sol_full, this->sol_full, this->sol_half);
    const real_t *err = this->convChecker->GetErrorNorms();

    // Calculate maximum error
    real_t maxErr = 0;
    len_t maxErri = 0;
    for (len_t i = 0; i < this->nontrivials.size(); i++) {
        const real_t scale = this->convChecker->GetErrorScale(this->nontrivials[i]);

        // scale = 0 indicates that |x| = 0, which could be ok
        if (scale != 0 && maxErr < err[i]/scale) {
            maxErr = err[i]/scale;
            maxErri = i;
        }
    }

    // Update time step
    const real_t S = 0.95;      // Safety factor
    const real_t alpha = (this->INCLUDE_STEP_STABILIZER ? 0.7 : 1.0);
    const real_t beta  = 0.4;

    real_t dt;
    if (this->constantStep) {      // "debug" mode
        dt = this->dt;
        converged = true;
    } else {
        dt = S*this->dt /pow(maxErr, alpha);

        if (this->INCLUDE_STEP_STABILIZER)
             dt *= pow(this->oldMaxErr, beta);
    }

    // Adjust for final time point
    real_t ctime = this->currentTime + this->oldDt;
    if (converged && ctime + dt > this->tMax) {
        dt = this->tMax - ctime;

        // Prevent round-off erors
        if (ctime+dt < this->tMax)
            dt *= (1 + std::numeric_limits<real_t>::epsilon());
    }

    if (this->verbose) {
        DREAM::IO::PrintInfo(
            "[TimeStepper] max error:  %.6e  (for unknown #" LEN_T_PRINTF_FMT ")\n"
            "[TimeStepper] step %s:  %.6e  ->  %.6e",
            maxErr, this->nontrivials[maxErri], ((dt < this->dt)?"DECREASED":"INCREASED"),
            this->dt, dt
        );
    }

    this->oldMaxErr = maxErr;
    this->dt = dt;

    return converged;
}

