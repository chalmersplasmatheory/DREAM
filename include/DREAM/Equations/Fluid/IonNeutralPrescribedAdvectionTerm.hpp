#ifndef _DREAM_EQUATION_ION_NEUTRAL_PRESCRIBED_ADVECTION_TERM_HPP
#define _DREAM_EQUATION_ION_NEUTRAL_PRESCRIBED_ADVECTION_TERM_HPP

#include "DREAM/Equations/Fluid/IonNeutralAdvectionTerm.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Interpolator1D.hpp"
#include "DREAM/MultiInterpolator1D.hpp"

namespace DREAM {
	class IonNeutralPrescribedAdvectionTerm : public IonNeutralAdvectionTerm {
	private:
	    len_t offset;
		MultiInterpolator1D *FrPrescribed;

	protected:
		virtual void SetPartialAdvectionTerm(len_t /*derivId*/, len_t /*nMultiples*/) override;
		virtual void SetCoeffs(const real_t) override;
		
	public:
		IonNeutralPrescribedAdvectionTerm(FVM::Grid *g, IonHandler *ihdl, const len_t iIon, bool allocCoefficients,
		 FVM::AdvectionInterpolationCoefficient::adv_interpolation intp, OptionConstants::adv_jacobian_mode jac_mode, 
            len_t id, real_t damping_factor, len_t offset, MultiInterpolator1D*);
		~IonNeutralPrescribedAdvectionTerm();
		
	};
}
#endif/*_DREAM_EQUATION_ION_NEUTRAL_PRESCRIBED_ADVECTION_TERM_HPP*/
