#ifndef _DREAM_EQUATION_ION_NEUTRAL_PRESCRIBED_DIFFUSION_TERM_HPP
#define _DREAM_EQUATION_ION_NEUTRAL_PRESCRIBED_DIFFUSION_TERM_HPP

#include "DREAM/Equations/Fluid/IonNeutralAdvectionDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Interpolator1D.hpp"
#include "DREAM/MultiInterpolator1D.hpp"

namespace DREAM {
	class IonNeutralPrescribedDiffusionTerm : public IonNeutralAdvectionDiffusionTerm<FVM::DiffusionTerm> {
	private:
	    len_t offset;
		MultiInterpolator1D *DrrPrescribed;

	protected:
		virtual void SetPartialDiffusionTerm(len_t /*derivId*/, len_t /*nMultiples*/) override;
		
	public:
		IonNeutralPrescribedDiffusionTerm(FVM::Grid *g, IonHandler *ihdl, const len_t iIon, bool allocCoefficients, len_t offset, MultiInterpolator1D*);
		~IonNeutralPrescribedDiffusionTerm();
		
		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
	};
}
#endif/*_DREAM_EQUATION_ION_NEUTRAL_PRESCRIBED_DIFFUSION_TERM_HPP*/
