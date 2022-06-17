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
	FVM::UnknownQuantityHandler *u
) : TimeStepper(u), tMax(tMax), dt0(dt0), dtMax(dtMax) {
	if (tMax < dt0)
		throw TimeStepperException("TimeStepperIonization: Maximum time must be greater than initial time step.");
	
	this->id_n_cold = u->GetUnknownID(OptionConstants::UQTY_N_COLD);
	this->nr = u->GetUnknown(this->id_n_cold)->GetGrid()->GetNr();

	this->timescales = new real_t[this->nr];
}

/**
 * Destructor.
 */
TimeStepperIonization::~TimeStepperIonization() {
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
	const real_t *np = n_cold->GetDataPrevious();

	real_t dt = this->currentTime - n_cold->GetPreviousTime();

	// Evaluate time scale as "1 / (|dn/dt| / n)
	for (len_t ir = 0; ir < this->nr; ir++)
		timescales[ir] = n[ir] * dt / std::abs(n[ir]-np[ir]);
	
	// Find shortest time scale
	real_t tscale = timescales[0];
	for (len_t ir = 0; ir < this->nr; ir++)
		if (timescales[ir] < tscale)
			tscale = timescales[ir];
	
	return tscale;
}

/**
 * Return the time for the next time step.
 */
real_t TimeStepperIonization::NextTime() {
	real_t timescale = GetIonizationTimeScale();
	
	if (this->tscale0 == 0) {
		this->tscale0 = timescale;
	}

	this->dt = (this->dt0 * timescale / this->tscale0);
	if (this->dtMax > 0 && this->dt > this->dtMax)
		this->dt = this->dtMax;

	return this->currentTime+this->dt;
}

/**
 * Validate the most recently taken time step. Currently, we
 * don't do this in this time stepper, but in the future we
 * could implement a check to ensure that the ionization time
 * scale is always respected and re-run if it is not.
 */
void TimeStepperIonization::ValidateStep() {
	this->currentTime += this->dt;
}

/**
 * Print the progress of the time stepper.
 */
void TimeStepperIonization::PrintProgress() {
	cout << "t = " << CurrentTime() << "s, dt = " << this->dt << "s" << endl;
}

