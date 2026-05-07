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
	FVM::UnknownQuantityHandler *unknowns, bool withIonJacobian
) : RechesterRosenbluthTransport(g, mgtype, dBB) {

	SetName("TrappingLimitedRRTransport");

	auto mg = g->GetMomentumGrid(0);
	if (mgtype != OptionConstants::MOMENTUMGRID_TYPE_PXI || mg->GetNp2() != 1)
		throw DREAMException("The 'TrappingLimitedRRTransport' operator is only applicable in isotropic (nxi = 1) mode.");

	this->nuD = cqh->GetNuD();

	this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
	this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
	this->id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

	AddUnknownForJacobian(unknowns, this->id_ncold);
	AddUnknownForJacobian(unknowns, this->id_Tcold);
	if (withIonJacobian)
		AddUnknownForJacobian(unknowns, this->id_ni);

	AllocateCache();
}


TrappingLimitedRRTransport::~TrappingLimitedRRTransport() {
	delete [] this->jacobianPrefactor;
}


void TrappingLimitedRRTransport::AllocateCache() {
	const len_t nr = this->grid->GetNr();

	delete [] this->jacobianPrefactor;

	this->np1_cache = this->grid->GetMomentumGrid(0)->GetNp1();
	this->jacobianPrefactor = new real_t[(nr+1) * this->np1_cache];
}


bool TrappingLimitedRRTransport::GridRebuilt() {
	this->FVM::DiffusionTerm::GridRebuilt();

	AllocateCache();
	return true;
}


/**
 * Rebuild the transport coefficient.
 */
void TrappingLimitedRRTransport::Rebuild(
	const real_t t, const real_t, FVM::UnknownQuantityHandler*
) {
	const real_t *dB_B = this->EvaluateDeltaBOverB(t);

	// TODO: Here I am assuming that all momentum grids are the same
	// at all radii and probably is always correct? 
	//
	const len_t nr = this->grid->GetNr();
	auto mg = this->grid->GetMomentumGrid(0);
	FVM::RadialGrid *rg = this->grid->GetRadialGrid();

	real_t qR0 = rg->GetR0();
	const real_t *p1 = mg->GetP1();
	const len_t np1 = this->np1_cache;
	real_t *const* nuD_f2 = this->nuD->GetValue_f2();
	const real_t *drf = this->deltaRadialFlux;

	// Major radius is set to 'inf' in cylindrical geometry.
	if (std::isinf(qR0))
		qR0 = 1.0;

	for (len_t ir = 0; ir < nr+1; ir++) {
		const real_t xiT = rg->GetXi0TrappedBoundary_fr(ir);
		const real_t fp = 0.5 * rg->GetFSA_B_f(ir) * (1 - xiT*xiT);
		const real_t ft = std::max<real_t>(0, 1 - fp);
		const real_t dBB2 = dB_B[ir] * dB_B[ir];
		const real_t w1 = (ir < nr) ? drf[ir] : 0.0;
		const real_t w2 = (ir > 0) ? 1.0 - drf[ir] : 0.0;
		real_t *prefactorRow = this->jacobianPrefactor + ir*np1;
		const bool trivial = (dBB2 == 0.0 || ft == 0.0 || xiT == 0.0);

		for (len_t i = 0; i < np1; i++) {
			const real_t p = p1[i];
			const real_t gamma = std::sqrt(1 + p*p);
			const real_t v = Constants::c * p / gamma;
			const real_t DrrIso = M_PI * qR0 * dBB2 * v * fp;

			if (trivial || v <= 0 || qR0 <= 0) {
				Drr(ir, i, 0) += DrrIso;
				prefactorRow[i] = 0.0;
				continue;
			}

			real_t nuD_f = 0.0;
			if (ir < nr)
				nuD_f += w1 * nuD_f2[ir][i];
			if (ir > 0)
				nuD_f += w2 * nuD_f2[ir-1][i];

			if (nuD_f <= 0) {
				Drr(ir, i, 0) += DrrIso;
				prefactorRow[i] = 0.0;
				continue;
			}

			const real_t tauStep = 2 * M_PI * qR0 / v;
			const real_t beta = xiT*xiT / (nuD_f * tauStep);
			const real_t h = 1 / (1 + ft * beta);

			Drr(ir, i, 0) += DrrIso * h;
			prefactorRow[i] = DrrIso * h * (1 - h) / nuD_f;
		}
	}
}


/**
 * Set Jacobian of the diffusion coefficients for this diffusion term.
 */
void TrappingLimitedRRTransport::SetPartialDiffusionTerm(len_t derivId, len_t nMultiples) {
	ResetDifferentiationCoefficients();

	if (
		derivId != this->id_ncold &&
		derivId != this->id_Tcold &&
		derivId != this->id_ni
	)
		return;

	// In isotropic PXi, p_f2 is cell-centered in momentum, so the j=0 P2 slice matches the p1[i] loop.
	const real_t *dNuD_f2 = this->nuD->GetUnknownPartialContribution(derivId, FVM::FLUXGRIDTYPE_P2);
	if (dNuD_f2 == nullptr)
		return;

	const len_t nr = this->grid->GetNr();
	auto mg = this->grid->GetMomentumGrid(0);
	const len_t np1 = this->np1_cache;
	const len_t rowStride = mg->GetNCells_f2();
	const len_t perMultiple = nr * rowStride;
	const real_t *drf = this->deltaRadialFlux;

	for (len_t n = 0; n < nMultiples; n++) {
		const real_t *dNuD_n = dNuD_f2 + n*perMultiple;
		for (len_t ir = 0; ir < nr+1; ir++) {
			const real_t w1 = (ir < nr) ? drf[ir] : 0.0;
			const real_t w2 = (ir > 0) ? 1.0 - drf[ir] : 0.0;
			const real_t *prefactorRow = this->jacobianPrefactor + ir*np1;

			for (len_t i = 0; i < np1; i++) {
				const real_t prefactor = prefactorRow[i];
				if (prefactor == 0.0)
					continue;

				real_t dNuD_f = 0.0;
				if (ir < nr)
					dNuD_f += w1 * dNuD_n[ir*rowStride + i];
				if (ir > 0)
					dNuD_f += w2 * dNuD_n[(ir-1)*rowStride + i];

				dDrr(ir, i, 0, n) = prefactor * dNuD_f;
			}
		}
	}
}
