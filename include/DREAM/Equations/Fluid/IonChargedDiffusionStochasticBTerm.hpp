#ifndef _DREAM_EQUATION_ION_CHARGED_DIFFUSION_STOCHASTIC_B_TERM_HPP
#define _DREAM_EQUATION_ION_CHARGED_DIFFUSION_STOCHASTIC_B_TERM_HPP

#include "DREAM/Equations/Fluid/IonChargedAdvectionDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Interpolator1D.hpp"
#include "DREAM/MultiInterpolator1D.hpp"

namespace DREAM {
	class IonChargedDiffusionStochasticBTerm : public IonChargedAdvectionDiffusionTerm<FVM::DiffusionTerm> {
	private:
		FVM::Interpolator1D *dBOverB;
		MultiInterpolator1D *DrrHat;
		
		real_t **dDrrdni;
		real_t **dDrrdWi;
		real_t **dDrrdNi;
		real_t **dDrrdTcold;
		
		len_t id_ni, id_Wi, id_Ni;
		
		void Allocate();
		void Deallocate();

	protected:
		virtual void SetCoeffs(const len_t Z0) override;
		virtual void SetCoeffsAllCS(const real_t dt) override;
		virtual void SetDiffCoeffsAllCS(const real_t t) override;
		virtual void SetPartialDiffusionTerm(len_t /*derivId*/, len_t /*nMultiples*/) override;
		
	public:
		IonChargedDiffusionStochasticBTerm(FVM::Grid *g, IonHandler *ihdl, bool allocCoefficients, const len_t iIon,
			FVM::Interpolator1D*, MultiInterpolator1D*, FVM::UnknownQuantityHandler *u);
		~IonChargedDiffusionStochasticBTerm();
	};
}
#endif/*_DREAM_EQUATION_ION_CHARGED_DIFFUSION_STOCHASTIC_B_TERM_HPP*/
