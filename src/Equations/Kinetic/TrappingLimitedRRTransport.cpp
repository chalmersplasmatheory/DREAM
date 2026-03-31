/**
 * Implementation of a Rechester-Rosenbluth diffusion coefficient which
 * is limited by the de-trapping time.
 */

#include "DREAM/Equations/Kinetic/TrappingLimitedRRTransport.hpp"


using namespace DREAM;


TrappingLimitedRRTransport::TrappingLimitedRRTransport(
	FVM::Grid *g, enum OptionConstants::momentumgrid_type mgtype,
	FVM::Interpolator1D *dBB, RunawayFluid *REFluid
) : RechesterRosenbluthTransport(g, mgtype, dBB), REFluid(REFluid) {

	SetName("TrappingLimitedRRTransport");

	auto mg = g->GetMomentumGrid(0);
	if (mg->GetNp2() != 1)
		throw DREAMException("The 'TrappingLimitedRRTransport' operator is only applicable in isotropic (nxi = 1) mode.");
}


TrappingLimitedRRTransport::~TrappingLimitedRRTransport() { }


/**
 * Rebuild the transport coefficient.
 */
void TrappingLimitedRRTransport::Rebuild() {
	const real_t *dB_B = this->deltaBOverB->Eval(t);

    // XXX Here we assume that all momentum grids are the same
    // at all radii...
    const len_t nr = this->grid->GetNr();
    auto mg = this->grid->GetMomentumGrid(0);

    real_t R0  = this->grid->GetRadialGrid()->GetR0();
    const real_t
        *p1 = mg->GetP1(),
        *p2 = mg->GetP2();
    const len_t
        np1 = mg->GetNp1(),
        np2 = mg->GetNp2();

    // Major radius is set to 'inf' in cylindrical geometry.
    // If so, we instead set it to 1 and effectively remove
    // it from the diffusion expression, allowing the user
    // to fully control the diffusion magnitude using the
    // dB/B parameter.
    if (isinf(R0))
        R0 = 1.0;

	PitchScatterFrequency *nuD = REFluid->GetNuD();
	real_t a = this->grid->GetRadialGrid->GetR_f()[nr];

    for (len_t ir = 0; ir < nr+1; ir++) {
        const real_t q = 1.0;   // TODO (safety factor)
        const real_t *BA_xi = grid->GetBA_xi_fr(ir);
        for (len_t i = 0; i < np1; i++) {
			// Evaluate D_loss
			real_t v = Constants::c * p1[i] / sqrt(1+p1[i]*p1[i]);
			real_t xiT = grid->GetRadialGrid()->GetXi0TrappedBoundary_fr(ir);
			real_t geometricFactor = 0.5 * grid->GetRadialGrid()->GetFSA_B_f(ir) * (1 - xiT*xiT);
			real_t Dloss = M_PI * q * R0 * dB_B[ir]*dB_B[ir] * v * geometricFactor;

			// Evaluate pitch-angle scattering coefficient
			Ddetrap = (a*a)/(xiT*xiT) * nuD->evaluateAtP(ir, p1[i], collSettings);

			Drr(ir,i,0) += std::min(Dloss, Ddetrap);
		}
	}
}


