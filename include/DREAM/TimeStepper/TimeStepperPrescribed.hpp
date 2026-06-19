#ifndef _DREAM_TIME_STEPPER_PRESCRIBED_HPP
#define _DREAM_TIME_STEPPER_PRESCRIBED_HPP

#include <iostream>
#include <vector>

#include "DREAM/TimeStepper/TimeStepper.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
/**
 * Time stepper which advances the simulation on a prescribed
 * array of time points:
 *
 *   t[0], t[1], ..., t[N]
 *
 * The number of time steps is Nt = N, and the step 'k' advances
 * from t[k-1] to t[k].
 *
 * NOTE: Input validation (monotonicity, size >= 2, etc.) is
 * intentionally not performed here and is expected to be done in
 * the SimulationGenerator/ConstructTimeStepper stage.
 */
class TimeStepperPrescribed : public TimeStepper {
   private:
    // Index of the *next* time point to advance to.
    // This follows the same convention as TimeStepperConstant:
    //   tIndex = 0   -> no steps taken yet
    //   after NextTime(): tIndex increments by 1
    len_t tIndex = 0;

    // Prescribed times
    std::vector<real_t> times;

    // Cached endpoints (for convenience/consistency with other steppers).
    real_t t0 = 0;
    real_t tMax = 0;

    // Total number of time steps: Nt = times.size()-1
    len_t Nt = 0;

    // Number of time steps to save to output (downsampling).
    // If 0: save every step.
    len_t nSaveSteps = 0, nextSaveStep_l = 0;
    real_t dSaveStep = 0, nextSaveStep = 0;

    void InitSaveSteps();

   public:
    /**
     * Constructor.
     *
     * ntimes:     Number of prescribed time points (>= 2).
     * times:      Array of prescribed time points (size = ntimes).
     * u:          Unknown quantity handler.
     * eqsys:      Equation system.
     * nSaveSteps: Number of internal steps to mark for saving (0 = save all).
     */
    TimeStepperPrescribed(
        const len_t ntimes, const real_t *times, FVM::UnknownQuantityHandler *u,
        EquationSystem *eqsys, const len_t nSaveSteps = 0
    );

    /**
     * Returns the time of the most recently completed time step.
     * If no steps have been completed, returns t[0].
     */
    virtual real_t CurrentTime() const override;

    /**
     * Returns 'true' if the time stepper has reached the final
     * prescribed time (or Python terminate hook requests exit).
     */
    virtual bool IsFinished() override;

    /**
     * Returns 'true' if this step should be saved to output.
     */
    virtual bool IsSaveStep() override;

    /**
     * Returns the maximum time for this simulation.
     */
    virtual real_t MaxTime() const override;

    /**
     * Advances internal step counter and returns the time of the next step.
     */
    virtual real_t NextTime() override;

    /**
     * Print current progress to stdout.
     */
    virtual void PrintProgress() override;

    /**
     * Validate the most recently taken time step.
     */
    virtual void ValidateStep() override;
};
}  // namespace DREAM

#endif /*_DREAM_TIME_STEPPER_PRESCRIBED_HPP*/
