/**
 * Implementation of the hyperresistive diffusion term with an adaptive
 * diffusion coefficient.
 */

#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/AdaptiveHyperresistiveDiffusionTerm.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
AdaptiveHyperresistiveDiffusionTerm::AdaptiveHyperresistiveDiffusionTerm(
	FVM::Grid *grid, FVM::UnknownQuantityHandler *uqh, IonHandler *ions,
	const real_t grad_j_tot_max, bool gradient_normalized,
	const real_t dBB0, bool localized
) : AdaptiveMHDLikeTransportTerm(grid, uqh, grad_j_tot_max, gradient_normalized, localized),
	HyperresistiveDiffusionTerm(grid, nullptr), ions(ions), dBB0(dBB0) {
	
	this->Lambda0 = new real_t[grid->GetNr()+1];
	this->Lambda = new real_t[grid->GetNr()+1];
	this->dLambda = new real_t[(grid->GetNr()+1) * ions->GetNzs()];

	this->id_n_i = uqh->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
}


/**
 * Destructor.
 */
AdaptiveHyperresistiveDiffusionTerm::~AdaptiveHyperresistiveDiffusionTerm() {
	delete [] this->Lambda;
	delete [] this->Lambda0;
	delete [] this->dLambda;
}


/**
 * Return the diffusion coefficient.
 *
 * t: Time to evaluate the coefficient at.
 */
const real_t *AdaptiveHyperresistiveDiffusionTerm::EvaluateLambda(const real_t t) {
	FVM::RadialGrid *rg = this->AdaptiveMHDLikeTransportTerm::grid->GetRadialGrid();
	const len_t nr = rg->GetNr();

	if (this->CheckTransportEnabled(t)) {
		const real_t R0 = rg->GetR0();

		if (isinf(R0)) {
			// In toroidal geometry, Lambda = inf
			for (len_t ir = 0; ir < nr+1; ir++)
				this->Lambda0[ir] = 100;
		} else {
			real_t qR0 = R0;
			// NOTE: Taking the average B is not entirely consistent, as
			// technically the entire coefficient Lambda should be flux-surface
			// averaged. However, taking dB/B as a constant, the coefficient
			// should evaluate to
			//
			//   <Lambda> = ... * <V_A> = ... * <B>.
			//
			for (len_t ir = 0; ir < nr+1; ir++) {
				real_t Bavg = rg->GetFSA_B(ir) * rg->GetBmin(ir);
				real_t V_A = Bavg / (sqrt(Constants::mu0) * this->ions->GetTotalIonMassDensity(ir));
				real_t a = rg->GetMinorRadius();
				// Total enclosed toroidal flux
				real_t psit = rg->GetToroidalFlux_f(nr);

				this->Lambda0[ir] =
					M_PI * Constants::mu0 * qR0 * R0 * V_A * (psit*psit)
					/ (144*a*a)
					* this->dBB0*this->dBB0;
			}
		}
	}
	
	for (len_t ir = 0; ir < nr+1; ir++)
		this->Lambda[ir] = Lambda0[ir] * this->mask[ir];

	// Set ion density derivative of Lambda
	const len_t nZs = this->ions->GetNzs();
	const real_t *n_i = uqh->GetUnknownData(this->id_n_i);
	for (len_t iZs = 0; iZs < nZs; iZs++) {
		for (len_t ir = 0; ir < nr+1; ir++) {
			this->dLambda[iZs*(nr+1) + ir] = -this->Lambda[ir] / (2*n_i[iZs*(nr+1) + ir]);
		}
	}
	
	return this->Lambda;
}


/**
 * Evaluate the partial derivative of Drr.
 */
void AdaptiveHyperresistiveDiffusionTerm::SetPartialDiffusionTerm(
	len_t derivId, len_t nMultiples
) {
	if (derivId != this->id_n_i)
		return;
	
	ResetDifferentiationCoefficients();

	for (len_t iZs = 0; iZs < nMultiples; iZs++)
		this->BuildCoefficient(this->dLambda + iZs*(nr+1), this->ddrr + iZs*(nr+1));
}


