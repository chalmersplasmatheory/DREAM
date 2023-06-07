/**
 * A variation in the magnitude of the magnetic field leads to an advective
 * process in pitch for electrons, due to the conservation of the magnetic
 * moment. This term approximately accounts for this fact and takes the
 * normalized rate-of-change of the magnetic field strength as input (assuming
 * that it is the toroidal magnetic field that is primarily being varied).
 */

#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/Kinetic/TimeVaryingBTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/fluxGridType.enum.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
TimeVaryingBTerm::TimeVaryingBTerm(
	FVM::Grid *g, FVM::Interpolator1D *dlnB0dt
) : FVM::AdvectionTerm(g), dlnB0dt(dlnB0dt) {
	
	SetName("TimeVaryingBTerm");

	BounceAverageForce();
}

/**
 * Destructor.
 */
TimeVaryingBTerm::~TimeVaryingBTerm() {
	for (len_t ir = 0; ir < this->nr; ir++)
		delete [] this->BA_Fxi[ir];
	delete [] this->BA_Fxi;
}

/**
 * Called whenever the grid is rebuilt (such as when the magnetic
 * field changes).
 */
bool TimeVaryingBTerm::GridRebuilt() {
	for (len_t ir = 0; ir < this->nr; ir++)
		delete [] this->BA_Fxi[ir];
	delete [] this->BA_Fxi;

	// Re-calculate bounce average of the compression force
	this->BounceAverageForce();

	return true;
}

/**
 * Bounce average the compressional force caused by the time-varying
 * magnetic field.
 */
void TimeVaryingBTerm::BounceAverageForce() {
	const len_t np1 = this->grid->GetMomentumGrid(0)->GetNp1();
	const len_t np2 = this->grid->GetMomentumGrid(0)->GetNp2();

	this->BA_Fxi = new real_t*[this->nr];
	this->BA_Fxi[0] = new real_t[this->nr*np1*(np2+1)];
	
	BAParams bap;
	for (len_t ir = 0; ir < this->nr; ir++) {
		FVM::MomentumGrid *mg = this->grid->GetMomentumGrid(ir);
		if (ir > 0)
			this->BA_Fxi[ir] = this->BA_Fxi[ir-1] + np1*(np2+1);

		for (len_t j_f = 0, idx = 0; j_f < np2+1; j_f++) {
			bap.xi0 = mg->GetP2_f(j_f);

			for (len_t i = 0; i < np1; i++, idx++) {
				BA_Fxi[ir][idx] = this->grid->CalculateBounceAverage(
					ir, i, j_f, FVM::FLUXGRIDTYPE_P2, BForce, &bap, nullptr
				);
			}
		}
	}
}

/**
 * Evaluate the (non-bounce averaged) compression force
 * caused by the time-varying magnetic field.
 */
real_t TimeVaryingBTerm::BForce(
	real_t xiOverXi0, real_t BOverBmin, real_t /*ROverR0*/,
	real_t /*NablaR2*/, void *par
) {
	BAParams *params = (BAParams*)par;
	
	real_t xi0 = params->xi0;
	real_t Fxi;
	if (xi0 != 0)
		//Fxi = BOverBmin*(1-xi0*xi0) / std::sqrt(1-BOverBmin*(1-xi0*xi0));
		Fxi = BOverBmin*(1-xi0*xi0) / (xiOverXi0*xi0);
	else
		Fxi = 0;
	
	return Fxi;
}

/**
 * Build the coefficients of this advection term.
 */
void TimeVaryingBTerm::Rebuild(
	const real_t t, const real_t,
	FVM::UnknownQuantityHandler*
) {
	real_t dBdt_B = this->dlnB0dt->Eval(t)[0];
	this->current_dlnB0dt = dBdt_B;

	const len_t np1 = this->grid->GetMomentumGrid(0)->GetNp1();
	const len_t np2 = this->grid->GetMomentumGrid(0)->GetNp2();

	for (len_t ir = 0; ir < this->nr; ir++) {
		for (len_t j_f = 0, idx = 0; j_f < np2+1; j_f++) {
			for (len_t i = 0; i < np1; i++, idx++) {
				F2(ir, i, j_f) -= BA_Fxi[ir][idx] * dBdt_B / 2;
			}
		}
	}
}

