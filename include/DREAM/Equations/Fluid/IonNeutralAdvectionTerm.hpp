#ifndef _DREAM_EQUATION_ION_NEUTRAL_ADVECTION_TERM_HPP
#define _DREAM_EQUATION_ION_NEUTRAL_ADVECTION_TERM_HPP

#include "DREAM/Equations/Fluid/IonNeutralAdvectionDiffusionTerm.hpp"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/AdvectionInterpolationCoefficient.hpp"

namespace DREAM {
	class IonNeutralAdvectionTerm : public IonNeutralAdvectionDiffusionTerm<FVM::AdvectionTerm> {
		
	public:
		IonNeutralAdvectionTerm(FVM::Grid *g, IonHandler *ihdl, const len_t iIon, bool allocCoefficients,
		FVM::AdvectionInterpolationCoefficient::adv_interpolation intp, OptionConstants::adv_jacobian_mode jac_mode, 
            len_t id, real_t damping_factor=1.0);
		~IonNeutralAdvectionTerm();
		
		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
		virtual void SetCoeffs(const real_t) = 0;
		
	};
}
#endif/*_DREAM_EQUATION_ION_NEUTRAL_ADVECTION_TERM_HPP*/
