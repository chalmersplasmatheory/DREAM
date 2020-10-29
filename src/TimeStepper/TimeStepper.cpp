/**
 * Common TimeStepper routines.
 */

#include "DREAM/TimeStepper/TimeStepper.hpp"


using namespace DREAM;


/**
 * Allocate memory for the solution vectors.
 *
 * size: Number of elements of full solution vector.
 */
void TimeStepper::AllocateSolutions(const len_t size) {
    if (this->sol_init != nullptr)
        this->DeallocateSolutions();

    this->sol_size = size;
    this->sol_init = new real_t[size];
}

/**
 * Deallocate memory for the solution vectors.
 */
void TimeStepper::DeallocateSolutions() {
    if (this->sol_init != nullptr) {
        delete [] this->sol_init;
        this->sol_init = nullptr;
    }
}

/**
 * Restore the initial solution.
 *
 * nSteps:   Number of time steps to roll back.
 * pushinit: If 'true', pushes the restored "initial" solution
 *           after rolling back.
 */
void TimeStepper::RestoreInitialSolution(const len_t nSteps, bool pushinit) {
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

