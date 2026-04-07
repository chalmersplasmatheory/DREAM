/**
 * Isotropic Rechester-Rosenbluth diffusion coefficient
 * limited by finite de-trapping time.
 */

#include <algorithm>
#include <cmath>
#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Kinetic/TrappingLimitedRRTransport.hpp"


using namespace DREAM;


TrappingLimitedRRTransport::TrappingLimitedRRTransport(
	FVM::Grid *g, enum OptionConstants::momentumgrid_type mgtype,
	FVM::Interpolator1D *dBB, CollisionQuantityHandler *cqh,
	FVM::UnknownQuantityHandler *unknowns
) : RechesterRosenbluthTransport(g, mgtype, dBB), unknowns(unknowns) {

	SetName("TrappingLimitedRRTransport");

	auto mg = g->GetMomentumGrid(0);
	if (mgtype != OptionConstants::MOMENTUMGRID_TYPE_PXI || mg->GetNp2() != 1)
		throw DREAMException("The 'TrappingLimitedRRTransport' operator is only applicable in isotropic (nxi = 1) mode.");

	this->nuD = cqh->GetNuD();
	this->collSettings = cqh->GetCollisionQuantitySettings();

	this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
	this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
	this->id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

	AddUnknownForJacobian(unknowns, this->id_ncold);
	AddUnknownForJacobian(unknowns, this->id_Tcold);
	AddUnknownForJacobian(unknowns, this->id_ni);
}


TrappingLimitedRRTransport::~TrappingLimitedRRTransport() { }


/**
 * Evaluate nu_D on the radial flux grid by linearly interpolating the two
 * neighboring radial cell values.
 */
real_t TrappingLimitedRRTransport::EvaluateNuDOnRadialFluxGrid(len_t ir, real_t p) const {
	const len_t nr = this->grid->GetNr();
	real_t nuD_f = 0;

	if (ir < nr)
		nuD_f += this->deltaRadialFlux[ir] * this->nuD->evaluateAtP(ir, p, this->collSettings);
	if (ir > 0)
		nuD_f += (1 - this->deltaRadialFlux[ir]) * this->nuD->evaluateAtP(ir-1, p, this->collSettings);

	return nuD_f;
}


/**
 * Evaluate the derivative of nu_D on the radial flux grid.
 */
real_t TrappingLimitedRRTransport::EvaluatePartialNuDOnRadialFluxGrid(
	len_t ir, real_t p, len_t derivId, len_t n
) const {
	const len_t nr = this->grid->GetNr();
	real_t dNuD_f = 0;

	if (ir < nr)
		dNuD_f += this->deltaRadialFlux[ir] * this->nuD->evaluatePartialAtP(ir, p, derivId, n, this->collSettings);
	if (ir > 0)
		dNuD_f += (1 - this->deltaRadialFlux[ir]) * this->nuD->evaluatePartialAtP(ir-1, p, derivId, n, this->collSettings);

	return dNuD_f;
}


/**
 * Rebuild the transport coefficient.
 */
void TrappingLimitedRRTransport::Rebuild(
	const real_t t, const real_t, FVM::UnknownQuantityHandler*
) {
	this->currentTime = t;
	const real_t *dB_B = this->EvaluateDeltaBOverB(t);

	// TODO: Here I am assuming that all momentum grids are the same
	// at all radii and probably is always correct? 
	//
	const len_t nr = this->grid->GetNr();
	auto mg = this->grid->GetMomentumGrid(0);
	FVM::RadialGrid *rg = this->grid->GetRadialGrid();

	real_t qR0 = rg->GetR0();
	const real_t *p1 = mg->GetP1();
	const len_t np1 = mg->GetNp1();

	// Major radius is set to 'inf' in cylindrical geometry.
	if (std::isinf(qR0))
		qR0 = 1.0;

	for (len_t ir = 0; ir < nr+1; ir++) {
		const real_t xiT = rg->GetXi0TrappedBoundary_fr(ir);
		const real_t fp = 0.5 * rg->GetFSA_B_f(ir) * (1 - xiT*xiT);
		const real_t ft = std::max<real_t>(0, 1 - fp);

		for (len_t i = 0; i < np1; i++) {
			const real_t p = p1[i];
			const real_t gamma = std::sqrt(1 + p*p);
			const real_t v = Constants::c * p / gamma;
			const real_t DrrIso = M_PI * qR0 * dB_B[ir]*dB_B[ir] * v * fp;

			if (DrrIso == 0 || ft == 0 || xiT == 0) {
				Drr(ir, i, 0) += DrrIso;
				continue;
			}

			const real_t nuD_f = EvaluateNuDOnRadialFluxGrid(ir, p);
			if (nuD_f <= 0 || v <= 0 || qR0 <= 0) {
				Drr(ir, i, 0) += DrrIso;
				continue;
			}

			const real_t tauStep = 2 * M_PI * qR0 / v;
			const real_t beta = xiT*xiT / (nuD_f * tauStep);
			const real_t h = 1 / (1 + ft * beta);

			Drr(ir, i, 0) += DrrIso * h;
		}
	}
}


/**
 * Set Jacobian of the diffusion coefficients for this diffusion term.
 */
void TrappingLimitedRRTransport::SetPartialDiffusionTerm(len_t derivId, len_t nMultiples) {
	ResetDifferentiationCoefficients();
	const real_t *dB_B = this->EvaluateDeltaBOverB(this->currentTime);

	if (
		derivId != this->id_ncold &&
		derivId != this->id_Tcold &&
		derivId != this->id_ni
	)
		return;

	const len_t nr = this->grid->GetNr();
	auto mg = this->grid->GetMomentumGrid(0);
	FVM::RadialGrid *rg = this->grid->GetRadialGrid();

	real_t qR0 = rg->GetR0();
	const real_t *p1 = mg->GetP1();
	const len_t np1 = mg->GetNp1();

	if (std::isinf(qR0))
		qR0 = 1.0;

	for (len_t n = 0; n < nMultiples; n++) {
		for (len_t ir = 0; ir < nr+1; ir++) {
			const real_t xiT = rg->GetXi0TrappedBoundary_fr(ir);
			const real_t fp = 0.5 * rg->GetFSA_B_f(ir) * (1 - xiT*xiT);
			const real_t ft = std::max<real_t>(0, 1 - fp);

			for (len_t i = 0; i < np1; i++) {
				const real_t p = p1[i];
				const real_t gamma = std::sqrt(1 + p*p);
				const real_t v = Constants::c * p / gamma;
				const real_t DrrIso = M_PI * qR0 * dB_B[ir]*dB_B[ir] * v * fp;

				if (DrrIso == 0 || ft == 0 || xiT == 0 || v <= 0 || qR0 <= 0)
					continue;

				const real_t nuD_f = EvaluateNuDOnRadialFluxGrid(ir, p);
				if (nuD_f <= 0)
					continue;

				const real_t tauStep = 2 * M_PI * qR0 / v;
				const real_t beta = xiT*xiT / (nuD_f * tauStep);
				const real_t h = 1 / (1 + ft * beta);
				const real_t dNuD_f = EvaluatePartialNuDOnRadialFluxGrid(ir, p, derivId, n);

				dDrr(ir, i, 0, n) = DrrIso * (1 - h) * dNuD_f / nuD_f;
			}
		}
	}
}
