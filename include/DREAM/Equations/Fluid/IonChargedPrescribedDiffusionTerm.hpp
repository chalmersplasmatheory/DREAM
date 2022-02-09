#ifndef _DREAM_EQUATION_ION_CHARGED_PRESCRIBED_DIFFUSION_TERM_HPP
#define _DREAM_EQUATION_ION_CHARGED_PRESCRIBED_DIFFUSION_TERM_HPP

#include "DREAM/Equations/Fluid/IonChargedAdvectionDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Interpolator1D.hpp"
#include "DREAM/MultiInterpolator1D.hpp"

namespace DREAM {
	class IonChargedPrescribedDiffusionTerm : public IonChargedAdvectionDiffusionTerm<FVM::DiffusionTerm> {
	private:
	    len_t offset;
		MultiInterpolator1D *DrrPrescribed;

	protected:
		virtual void SetCoeffs(const len_t Z0) override;
		virtual void SetCoeffsAllCS(const real_t dt) override;
		virtual void SetDiffCoeffsAllCS(const real_t t) override;
		virtual void SetPartialDiffusionTerm(len_t /*derivId*/, len_t /*nMultiples*/) override;
		
	public:
		IonChargedPrescribedDiffusionTerm(FVM::Grid *g, IonHandler *ihdl, const len_t iIon, bool allocCoefficients, len_t offset, MultiInterpolator1D*);
		~IonChargedPrescribedDiffusionTerm();
	};
}
#endif/*_DREAM_EQUATION_ION_CHARGED_PRESCRIBED_DIFFUSION_TERM_HPP*/
