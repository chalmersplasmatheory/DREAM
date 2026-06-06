/**
 * Implementation of the hyperresistive diffusion of poloidal flux.
 * Based on the section on hyperresistivity in DREAM/doc/notes/theory
 */

#include "DREAM/Constants.hpp"
#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/Fluid/HyperresistiveDiffusionTerm.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
HyperresistiveDiffusionTerm::HyperresistiveDiffusionTerm(
    FVM::Grid *g, IonHandler *ions, FVM::Interpolator1D *Lambda, FVM::Interpolator1D *dBB
) : FVM::DiffusionTerm(g), ions(ions), Lambda(Lambda), dBB(dBB) {
    
    SetName("HyperresistiveDiffusionTerm");

	if (this->dBB != nullptr)
		this->lambda_buf = new real_t[g->GetNr()+1];
}


/**
 * Destructor.
 */
HyperresistiveDiffusionTerm::~HyperresistiveDiffusionTerm() {
	if (this->Lambda != nullptr)
		delete this->Lambda;
	if (this->dBB != nullptr)
		delete this->dBB;
	if (this->lambda_buf != nullptr)
		delete [] this->lambda_buf;
}


/**
 * Get the current value for the hyperresistive diffusion coefficient.
 */
const real_t *HyperresistiveDiffusionTerm::GetLambda() const {
	if (this->Lambda != nullptr)
		return this->Lambda->GetBuffer();
	else if (this->dBB != nullptr)
		return this->lambda_buf;
	else
		throw DREAMException("HyperresistiveDiffusionTerm: either Lambda or dBB must be provided.");
}


/**
 * Evaluate the hyperresistive diffusion coefficient.
 */
const real_t *HyperresistiveDiffusionTerm::EvaluateLambda(const real_t t) {
	if (this->Lambda != nullptr)
		return this->Lambda->Eval(t);
	else if (this->dBB != nullptr) {
		// Evaluate Lambda from dB/B
		FVM::RadialGrid *rGrid = grid->GetRadialGrid(); 
		const len_t nr = rGrid->GetNr();
		const real_t R0 = rGrid->GetR0();

		const real_t *dBB0 = this->dBB->Eval(t);
		const real_t qR0 = R0;

		// (skip ir=0 since psi_t=0 there, and to avoid 1/psitPrime = 1/0)
		for (len_t ir = 0; ir < nr+1; ir++) {
			real_t Bavg = rGrid->GetFSA_B_f(ir) * rGrid->GetBmin_f(ir);

			real_t rho;
			if (ir < nr)
				rho = this->ions->GetTotalIonMassDensity(ir);
			else
				rho = this->ions->GetTotalIonMassDensity(nr-1);

			real_t V_A = Bavg / std::sqrt(Constants::mu0 * rho);
			real_t a = rGrid->GetMinorRadius();
			real_t psit = rGrid->GetToroidalFlux_f(nr);

			real_t pf = M_PI*Constants::mu0 * qR0 * R0;

			this->lambda_buf[ir] =
				pf * V_A * (psit*psit) / (288*a*a) * dBB0[ir]*dBB0[ir];
		}

		return this->lambda_buf;
	} else
		throw DREAMException("HyperresistiveDiffusionTerm: either Lambda or dBB must be provided.");
}


/**
 * Build the coefficients of this diffusion term.
 */
void HyperresistiveDiffusionTerm::Rebuild(
	const real_t t, const real_t, FVM::UnknownQuantityHandler*
) {
    const real_t *Lmbd  = this->EvaluateLambda(t);
	this->BuildCoefficient(Lmbd, this->drr);
}

/**
 * Set diffusion coefficient, or derivative of, diffusion coefficient.
 */
void HyperresistiveDiffusionTerm::BuildCoefficient(
	const real_t *coeff, real_t **diffusion_coeff
) {
    FVM::RadialGrid *rGrid = grid->GetRadialGrid(); 

    // XXX: here we assume that all radii have the same momentum grids
    const len_t np1 = n1[0], np2 = n2[0];

    // (skip ir=0 since psi_t=0 there, and to avoid 1/psitPrime = 1/0)
    for (len_t ir = 1; ir < nr+1; ir++) {
        real_t Bmin = rGrid->GetBmin_f(ir);
        real_t BdotPhi = rGrid->GetBTorG_f(ir)*rGrid->GetFSA_1OverR2_f(ir);
        real_t VpVol = rGrid->GetVpVol_f(ir); 

        real_t psitPrime = VpVol*BdotPhi / (2*M_PI);

        // The entire d psi/dt equation is multiplied by 2*pi*psi_t'/VpVol,
        // so we get an extra factor of 2*pi here, and only one factor
        // of psi_t' (the factor 1/VpVol from the normalisation is included 
        // via the DiffusionTerm class). We also add a factor of 1/VpVol 
        // to the diffusion coefficient to cancel the factor VpVol inside 
        // the first radial derivative added by the DiffusionTerm. 
        //
        // Also, we divide by 'Bmin' since this operator is applied to
        // 'j_tot / (B/Bmin)'.
        real_t drr = 
            2*M_PI*rGrid->GetToroidalFlux_f(ir)*coeff[ir] / (VpVol*psitPrime*Bmin);

        for (len_t j = 0; j < np2; j++) 
            for (len_t i = 0; i < np1; i++) 
				diffusion_coeff[ir][j*np1 + i] += drr;
    }
}

