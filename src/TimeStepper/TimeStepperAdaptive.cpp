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
 * tMax:        Final simulation time.
 * dt0:         Initial time step.
 * uqh:         UnknownQuantityHandler of solver.
 * nontrivials: List of non-trivial unknowns.
 * reltol:      Default relative tolerance.
 * checkEvery:  Number of time steps to take _without_ doing a convergence
 *              check after each check (i.e. 0 => check _every_ time step,
 *              1 => check every other etc.)
 */
TimeStepperAdaptive::TimeStepperAdaptive(
    const real_t tMax, const real_t dt0, FVM::UnknownQuantityHandler *uqh,
    vector<len_t>& nontrivials, const real_t reltol, int_t checkEvery,
    bool verbose
) : TimeStepper(uqh), tMax(tMax), dt(dt0), nontrivials(nontrivials),
  checkEvery(checkEvery), verbose(verbose) {
    
    this->stepsSinceCheck = checkEvery;
    this->convChecker = new ConvergenceChecker(uqh, nontrivials, reltol);

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
    DeallocateSolutionFull();
    DeallocateSolutions();
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
    
    switch (stg) {
        case STAGE_NORMAL:
            if (ShouldCheckError()) {
                this->stepsSinceCheck = 0;
                stg = STAGE_FIRST_HALF;
            } else {
                this->stepsSinceCheck++;
                this->currentTime += this->dt;
            }
            break;

        case STAGE_FIRST_HALF:
            stg = STAGE_SECOND_HALF;
            this->initTime = this->currentTime;
            this->currentTime += 0.5*this->dt;
            break;

        case STAGE_SECOND_HALF:
            stg = STAGE_FULL;
            // Get solution due to two half steps
            CopySolution(&this->sol_half);
            RestoreInitialSolution();
            this->currentTime = this->initTime;
            break;

        case STAGE_FULL: {
            // This means that the previous step was the full step,
            // so that we should now evaluate how successful the time
            // stepping has been. We must also update the current time
            // (and make sure it's done before modifying the time step dt)
            real_t oldDt = this->dt;
            if (UpdateStep()) {
                stg = (checkEvery > 0) ? STAGE_NORMAL : STAGE_FIRST_HALF;
                this->currentTime += oldDt;
            } else {  // Time-stepping failed => redo step
                stg = STAGE_FIRST_HALF;
                RestoreInitialSolution();
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
 * Allocate memory for storing the initial solution (which contains
 * data for _all_ unknowns).
 */
void TimeStepperAdaptive::AllocateSolutionFull(const len_t size) {
    if (this->sol_init != nullptr)
        DeallocateSolutionFull();

    this->sol_init_size = size;
    this->sol_init = new real_t[size];
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
    this->sol_half = new real_t[size];
    this->sol_full = new real_t[size];
}

/**
 * Deallocate memory for the initial solution vector.
 */
void TimeStepperAdaptive::DeallocateSolutionFull() {
    if (this->sol_init != nullptr)
        delete [] this->sol_init;
}

/**
 * Deallocate memory for the solution vectors.
 */
void TimeStepperAdaptive::DeallocateSolutions() {
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
 * Copy full solution vector to temporary storage.
 */
void TimeStepperAdaptive::CopySolutionFull(real_t **sol) {
    const len_t FSIZE = this->unknowns->GetLongVectorSizeAll();
    if (sol_init_size != FSIZE)
        AllocateSolutionFull(FSIZE);

    this->unknowns->GetLongVectorAll(*sol);
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
 * Returns 'true' if the simulation has reached the final
 * time. Returns 'false' otherwise.
 */
bool TimeStepperAdaptive::IsFinished() {
    return 
        ((this->currentStage == STAGE_NORMAL || this->currentStage == STAGE_FIRST_HALF)
        && this->currentTime >= this->tMax);
}

/**
 * Returns 'true' if the result of the current time step
 * should be saved to the output.
 */
bool TimeStepperAdaptive::IsSaveStep() {
    return (this->currentStage == STAGE_NORMAL || this->currentStage == STAGE_FULL);
}

/**
 * Calculate and return the time of the next step to take.
 * The asbolute time value is returned and the caller will
 * have to turn that into a time step themselves.
 */
real_t TimeStepperAdaptive::NextTime() {
    real_t newTime = this->initTime;

    ts_stage stg = AdvanceStage();

    switch (stg) {
        case STAGE_NORMAL:
            newTime += this->dt;
            break;

        case STAGE_FIRST_HALF:
            this->CopySolutionFull(&this->sol_init);
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
    if (this->currentStage != STAGE_NORMAL)
        return;

    const len_t PERC_FMT_PREC = 2;      // Precision (after decimal point) in percentage
    //                          100 . XX            %
    const len_t PERC_FMT_LENGTH = 3+1+PERC_FMT_PREC+1;
    const len_t EDGE_LENGTH = 1;
    const len_t PROG_LENGTH = PROGRESSBAR_LENGTH-2*EDGE_LENGTH - PERC_FMT_LENGTH - 1;

    cout << "\r[";
    real_t perc     = CurrentTime()/this->tMax;
    len_t threshold = static_cast<len_t>(perc * PROG_LENGTH);
    
    for (len_t i = 0; i < PROG_LENGTH; i++) {
        if (i < threshold)
            cout << '=';
        else
            cout << '-';
    }

    cout << "] ";
    printf("%*.*f%%", int(4+PERC_FMT_PREC), int(PERC_FMT_PREC), perc*100.0);

    // Ensure that output is written right away (otherwise it may
    // not be written until the end-of-line character is written)
    cout << flush;
}

/**
 * Restore the initial solution.
 */
void TimeStepperAdaptive::RestoreInitialSolution() {
    this->unknowns->SetFromLongVectorAll(this->sol_init, true);
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
    CopySolution(&this->sol_full);

    bool converged = this->convChecker->IsConverged(this->sol_full, this->sol_full, this->sol_half);
    const real_t *err = this->convChecker->GetErrorNorms();

    // Calculate maximum error
    real_t maxErr = 0;
    for (len_t i = 0; i < this->nontrivials.size(); i++) {
        const real_t scale = this->convChecker->GetErrorScale(this->nontrivials[i]);

        // scale = 0 indicates that |x| = 0, which could be ok
        if (scale != 0 && maxErr < err[i]/scale)
            maxErr = err[i]/scale;
    }

    // Update time step
    const real_t S = 0.95;      // Safety factor
    const real_t alpha = 0.7;
    const real_t beta  = 0.4;

    real_t dt = S*this->dt*pow(this->oldMaxErr, beta)/pow(maxErr, alpha);

    if (this->verbose) {
        DREAM::IO::PrintInfo(
            "[TimeStepper] error:    %.6e\n"
            "[TimeStepper] step %s:  %.6e  ->  %.6e",
            maxErr, (maxErr>1?"DECREASED":"INCREASED"),
            this->dt, dt
        );
    }

    this->oldMaxErr = maxErr;
    this->dt = dt;

    return converged;
}

