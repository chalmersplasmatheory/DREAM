/**
 * Construct a time stepper object.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/TimeStepper/TimeStepper.hpp"
#include "DREAM/TimeStepper/TimeStepperAdaptive.hpp"
#include "DREAM/TimeStepper/TimeStepperConstant.hpp"
#include "DREAM/TimeStepper/TimeStepperIonization.hpp"
#include "DREAM/TimeStepper/TimeStepperPrescribed.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM;
using namespace std;


#define MODULENAME "timestep"


/**
 * Define options for the time stepper.
 *
 * s: Settings object to define settings in.
 */
void SimulationGenerator::DefineOptions_TimeStepper(Settings *s) {
	s->DefineSetting(MODULENAME "/automaticstep", "Step length for the automatic determination of the time step in the ionization time stepper.", (real_t)1e-12);
    s->DefineSetting(MODULENAME "/checkevery", "Check the error every N'th step (0 = check error after _every_ time step)", (int_t)0);
    s->DefineSetting(MODULENAME "/constantstep", "Override the adaptive stepper and force a constant time step (DEBUG OPTION)", (bool)false);
    s->DefineSetting(MODULENAME "/dt", "Length of each time step", (real_t)0.0);
	s->DefineSetting(MODULENAME "/dtmax", "Maximum allowed time step for the adaptive ionization time stepper.", (real_t)0);
	s->DefineSetting(MODULENAME "/minsavedt", "Minimum time required to elapse between time steps to save in adaptive ioniz time stepper.", (real_t)0);
    s->DefineSetting(MODULENAME "/nsavesteps", "Number of time steps to save to output (downsampling)", (int_t)0);
    s->DefineSetting(MODULENAME "/nt", "Number of time steps to take", (int_t)0);
	s->DefineSetting(MODULENAME "/safetyfactor", "Safety factor to use when automatically determining the baseline timestep for the adaptive ionization time stepper.", (real_t)50);
    s->DefineSetting(MODULENAME "/tmax", "Maximum simulation time", (real_t)0.0);
	s->DefineSetting(MODULENAME "/alpha", "Scaling factor for the ionization time. If not zero, timesteps are updated as: dt = min(dtmax, alpha * t_ioniz)", (real_t)0.0);
    s->DefineSetting(MODULENAME "/type", "Time step generator type", (int_t)OptionConstants::TIMESTEPPER_TYPE_CONSTANT);
    s->DefineSetting(MODULENAME "/verbose", "If true, generates excessive output", (bool)false);
    s->DefineSetting(MODULENAME "/times", "Prescribed simulation time points (1D array). The simulation advances from times[k-1] to times[k].", (len_t)0, (const real_t*)nullptr);
#ifdef DREAM_IS_PYTHON_LIBRARY
    s->DefineSetting(MODULENAME "/terminatefunc", "Python function used to determine when to terminate time stepping", (void*)nullptr);
#endif

    // Tolerance settings for adaptive time stepper
    DefineToleranceSettings(MODULENAME, s);
}

/**
 * Construct a TimeStepper object according to the settings.
 *
 * eqsys: Equation system object to assign TimeStepper object to.
 * s:     Settings object specifying how to construct the TimeStepper object.
 */
void SimulationGenerator::ConstructTimeStepper(EquationSystem *eqsys, Settings *s) {
    enum OptionConstants::timestepper_type type = (enum OptionConstants::timestepper_type)s->GetInteger(MODULENAME "/type");

    FVM::UnknownQuantityHandler *u = eqsys->GetUnknownHandler();
    vector<len_t> *nontrivials = eqsys->GetNonTrivialUnknowns();
    TimeStepper *ts;
    switch (type) {
        case OptionConstants::TIMESTEPPER_TYPE_CONSTANT:
            ts = ConstructTimeStepper_constant(s, u, eqsys);
            break;

        case OptionConstants::TIMESTEPPER_TYPE_ADAPTIVE:
            ts = ConstructTimeStepper_adaptive(s, u, eqsys, nontrivials);
            break;

		case OptionConstants::TIMESTEPPER_TYPE_IONIZATION:
			ts = ConstructTimeStepper_ionization(s, u, eqsys);
			break;

        case OptionConstants::TIMESTEPPER_TYPE_PRESCRIBED:
            ts = ConstructTimeStepper_prescribed(s, u, eqsys);
            break;

        default:
            throw SettingsException(
                "Unrecognized time stepper type: %d.", type
            );
    }

#ifdef DREAM_IS_PYTHON_LIBRARY
    void *terminatefunc = s->GetAddress(MODULENAME "/terminatefunc");
    if (terminatefunc != nullptr)
        ts->SetPythonTerminateFunc(terminatefunc);
#endif

    eqsys->SetTimeStepper(ts);
}


/**
 * Construct a TimeStepperConstant object according to the
 * provided settings.
 *
 * s: Settings object specifying how to construct the
 *    TimeStepperConstant object.
 */
TimeStepperConstant *SimulationGenerator::ConstructTimeStepper_constant(
    Settings *s, FVM::UnknownQuantityHandler *u,
    EquationSystem *eqsys
) {
    real_t tmax = s->GetReal(MODULENAME "/tmax");
    real_t dt   = s->GetReal(MODULENAME "/dt", false);
    int_t nt    = s->GetInteger(MODULENAME "/nt", false);
    int_t nSaveSteps = s->GetInteger(MODULENAME "/nsavesteps");

    bool dtset  = (dt > 0);
    bool ntset  = (nt > 0);

    if (dtset && ntset)
        throw SettingsException(
            "TimeStepper constant: "
            "Ambiguous time step specified. Only one of 'dt' and 'nt' "
            "may be set for the time stepper."
        );
    else if (!dtset && !ntset)
        throw SettingsException(
            "TimeStepper constant: "
            "No time step specified. Exactly one of 'dt' and 'nt' "
            "must be set for the time stepper."
        );

    // Generate object
    if (dtset) {
        s->MarkUsed(MODULENAME "/dt");
        return new TimeStepperConstant(tmax, dt, u, eqsys, nSaveSteps);
    } else {
        s->MarkUsed(MODULENAME "/nt");
        return new TimeStepperConstant(tmax, (len_t)nt, u, eqsys, nSaveSteps);
    }
}

/**
 * Construct a TimeStepperAdaptive object according to the
 * provided settings.
 *
 * s: Settings object specifying how to construct the
 *    TimeStepperAdaptive object.
 */
TimeStepperAdaptive *SimulationGenerator::ConstructTimeStepper_adaptive(
    Settings *s, FVM::UnknownQuantityHandler *u,
	EquationSystem *eqsys, vector<len_t> *nontrivials
) {
    int_t checkevery = s->GetInteger(MODULENAME "/checkevery");
    real_t tmax = s->GetReal(MODULENAME "/tmax");
    real_t dt = s->GetReal(MODULENAME "/dt");
    bool verbose = s->GetBool(MODULENAME "/verbose");
    bool conststep = s->GetBool(MODULENAME "/constantstep");

	vector<UnknownQuantityEquation*> *eqns = eqsys->GetEquations();

    if (dt == 0)
        dt = 1;

    ConvergenceChecker *cc = LoadToleranceSettings(
        MODULENAME, s, eqns, u, *nontrivials
    );

    return new TimeStepperAdaptive(tmax, dt, u, eqsys, *nontrivials, eqns, cc, checkevery, verbose, conststep);
}

/**
 * Construct a TimeStepperIonization object according to the
 * provided settings.
 *
 * s: Settings object specifying how to construct the
 *    TimeStepperIonization object.
 */
TimeStepperIonization *SimulationGenerator::ConstructTimeStepper_ionization(
	Settings *s, FVM::UnknownQuantityHandler *u, EquationSystem *eqsys
) {
	real_t automaticstep = s->GetReal(MODULENAME "/automaticstep");
	real_t dt = s->GetReal(MODULENAME "/dt");
	real_t dtmax = s->GetReal(MODULENAME "/dtmax");
	real_t minSaveDt = s->GetReal(MODULENAME "/minsavedt");
	real_t safetyfactor = s->GetReal(MODULENAME "/safetyfactor");
	real_t tmax = s->GetReal(MODULENAME "/tmax");
	real_t alpha = s->GetReal(MODULENAME "/alpha");

	if (dt < 0)
		throw SettingsException("TimeStepper ionization: Initial time step 'dt0' must be non-negative.");

	return new TimeStepperIonization(tmax, dt, dtmax, u, eqsys, automaticstep, safetyfactor, minSaveDt, alpha);
}

/**
 * Construct a TimeStepperPrescribed object according to the
 * provided settings.
 *
 * s: Settings object specifying how to construct the
 *    TimeStepperPrescribed object.
 */
TimeStepperPrescribed *SimulationGenerator::ConstructTimeStepper_prescribed(
    Settings *s, FVM::UnknownQuantityHandler *u,
    EquationSystem *eqsys
) {
    // tmax is still part of the standard interface. For prescribed stepping,
    // we require that (if set) it matches the last time point.
    real_t tmax = s->GetReal(MODULENAME "/tmax");

    // Disallow mixing with dt/nt for clarity.
    real_t dt = s->GetReal(MODULENAME "/dt", false);
    int_t  nt = s->GetInteger(MODULENAME "/nt", false);
    if (dt > 0 || nt > 0)
        throw SettingsException(
            "TimeStepper prescribed: "
            "Ambiguous time step specified. 'dt'/'nt' cannot be used together "
            "with the prescribed time grid 'times'."
        );

    int_t nSaveSteps_i = s->GetInteger(MODULENAME "/nsavesteps");
    if (nSaveSteps_i < 0)
        throw SettingsException(
            "TimeStepper prescribed: "
            "Invalid negative value assigned to 'nsavesteps': %d.", nSaveSteps_i
        );
    len_t nSaveSteps = (len_t)nSaveSteps_i;

    // Load prescribed time array (1D).
    len_t dims[1] = {0};
    const real_t *times = s->GetRealArray(MODULENAME "/times", 1, dims);
    const len_t ntimes = dims[0];

    if (ntimes < 2)
        throw SettingsException(
            "TimeStepper prescribed: "
            "Invalid time grid. 'times' must contain at least two time points."
        );

    // Validate time grid (explicit, no cleverness).
    if (times[0] != 0)
        throw SettingsException(
            "TimeStepper prescribed: "
            "Invalid time grid. The first time point must be 0."
        );

    for (len_t i = 1; i < ntimes; i++) {
        if (times[i] <= times[i-1])
            throw SettingsException(
                "TimeStepper prescribed: "
                "Invalid time grid. 'times' must be strictly increasing."
            );
    }

    const real_t tmax_times = times[ntimes-1];
    if (tmax <= 0) {
        // If user did not set tmax, accept and implicitly use the last time point.
        // (We still keep 'tmax' in settings for consistency with other steppers.)
        tmax = tmax_times;
    } else if (tmax != tmax_times) {
        throw SettingsException(
            "TimeStepper prescribed: "
            "Inconsistent final time. 'tmax' (= %.16e) must match times[end] (= %.16e).",
            tmax, tmax_times
        );
    }

    return new TimeStepperPrescribed(ntimes, times, u, eqsys, nSaveSteps);
}
