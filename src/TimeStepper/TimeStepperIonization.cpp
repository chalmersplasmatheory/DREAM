/**
 * Implementation of an adaptive time stepper which gradually
 * rescales the time step to ensure that the ionization time
 * scale is (just barely) respected.
 */

#include <iostream>
#include "DREAM/IO.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/TimeStepper/TimeStepperIonization.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
TimeStepperIonization::TimeStepperIonization(
	const real_t tMax, const real_t dt0, const real_t dtMax,
	FVM::UnknownQuantityHandler *u, EquationSystem *eqsys,
	const real_t automaticstep, const real_t safetyfactor,
	const real_t minSaveDt
) : TimeStepper(u, eqsys), tMax(tMax), dt0(dt0), dtMax(dtMax),
	minSaveDt(minSaveDt), automaticstep(automaticstep),
	safetyfactor(safetyfactor)
{
	if (tMax < dt0)
		throw TimeStepperException("TimeStepperIonization: Maximum time must be greater than initial time step.");
	
	this->id_n_cold = u->GetUnknownID(OptionConstants::UQTY_N_COLD);
	this->nr = u->GetUnknown(this->id_n_cold)->GetGrid()->GetNr();

	this->timescales = new real_t[this->nr];
	this->ncold = new real_t[this->nr];
}

/**
 * Destructor.
 */
TimeStepperIonization::~TimeStepperIonization() {
	delete [] this->ncold;
	delete [] this->timescales;
}


/**
 * Evaluates the current ionization time scale. If the time
 * scale cannot be evaluated (i.e. if n_cold is only known in
 * a single time point), zero is returned.
 */
real_t TimeStepperIonization::GetIonizationTimeScale() {
	if (this->currentTime == 0)
		return 0;

	FVM::UnknownQuantity *n_cold = this->unknowns->GetUnknown(this->id_n_cold);
	const real_t *n  = n_cold->GetData();

	// Evaluate time scale as "1 / (|dn/dt| / n)"
	for (len_t ir = 0; ir < this->nr; ir++) {
		if (n[ir] != this->ncold[ir])
			timescales[ir] = n[ir] * this->dt / std::abs(n[ir]-this->ncold[ir]);
		else
			timescales[ir] = std::numeric_limits<real_t>::infinity();
	}
	
	// Find shortest time scale
	real_t tscale = timescales[0];
	for (len_t ir = 0; ir < this->nr; ir++)
		if (timescales[ir] < tscale)
			tscale = timescales[ir];
	
	if (isinf(tscale))
		return 0;
	else
		return tscale;
}

/**
 * Return the time for the next time step.
 */
real_t TimeStepperIonization::NextTime() {
	real_t timescale = GetIonizationTimeScale();

	if (timescale == 0) {
		if (this->dt0 == 0)
			this->dt = this->automaticstep;
		else
			this->dt = this->dt0;
	} else {
		if (this->tscale0 == 0)
			this->tscale0 = timescale;

		// If the baseline time step should be set automatically...
		if (this->dt0 == 0) {
			this->dt0 = timescale / this->safetyfactor;
			DREAM::IO::PrintInfo("Setting baseline timestep: %.7e\n", this->dt0);
		}

		this->dt = (this->dt0 * timescale / this->tscale0);
		if (this->dtMax > 0 && this->dt > this->dtMax)
			this->dt = this->dtMax;
	}

	if ((this->currentTime+this->dt) > this->tMax)
		return this->tMax;
	else
		return this->currentTime+this->dt;
}

/**
 * Returns the maximum time for this simulation.
 */
real_t TimeStepperIonization::MaxTime() const {
    return this->tMax;
}

/**
 * Returns 'true' if the current time step is to be saved to the output file.
 */
bool TimeStepperIonization::IsSaveStep() {
	if (this->minSaveDt == 0)
		return true;
	else if (
		(this->lastSaveTime+this->minSaveDt) < this->currentTime ||
		this->tMax <= this->currentTime
	) {
		this->lastSaveTime = this->currentTime;
		return true;
	} else
		return false;
}

/**
 * Validate the most recently taken time step. Currently, we
 * don't do this in this time stepper, but in the future we
 * could implement a check to ensure that the ionization time
 * scale is always respected and re-run if it is not.
 */
void TimeStepperIonization::ValidateStep() {
	if ((this->currentTime+this->dt) > this->tMax)
		this->currentTime = this->tMax;
	else
		this->currentTime += this->dt;
	
	this->currentStep++;

	// Copy current value of n_cold
	const real_t *n = this->unknowns->GetUnknownDataPrevious(this->id_n_cold);
	for (len_t ir = 0; ir < this->nr; ir++)
		this->ncold[ir] = n[ir];
}

/**
 * Print the progress of the time stepper.
 */
void TimeStepperIonization::PrintProgress() {
	const len_t PERC_FMT_PREC = 2; // Precision (after decimal point) in percent
	//                          100 . XX            %
	const len_t PERC_FMT_LENGTH = 3+1+PERC_FMT_PREC+1;
	const len_t EDGE_LENGTH = 1;
	const len_t PROG_LENGTH =
		PROGRESSBAR_LENGTH - 2*EDGE_LENGTH - PERC_FMT_LENGTH - 1;
	
	cout << "\r[";
	real_t perc     = CurrentTime()/this->tMax;
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

	cout << flush;
}

