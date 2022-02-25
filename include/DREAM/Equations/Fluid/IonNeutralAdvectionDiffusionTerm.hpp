#ifndef _DREAM_EQUATION_ION_NEUTRAL_ADVECTION_DIFFUSION_TERM_HPP
#define _DREAM_EQUATION_ION_NEUTRAL_ADVECTION_DIFFUSION_TERM_HPP

#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"

namespace DREAM {
	template<class T>
    class IonNeutralAdvectionDiffusionTerm : public IonEquationTerm<T> {
		
	public:
		IonNeutralAdvectionDiffusionTerm(FVM::Grid *g, IonHandler *ihdl, const len_t iIon, bool allocCoefficients);
		virtual ~IonNeutralAdvectionDiffusionTerm();
		
		virtual bool SetCSJacobianBlock(
		    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *nions,
		    const len_t iIon, const len_t Z0, const len_t rOffset
		) override;
		virtual void SetCSMatrixElements(
		    FVM::Matrix *mat, real_t *rhs, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
		) override;
		virtual void SetCSVectorElements(
		    real_t *vec, const real_t *nions, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
		) override;
	};
	
	#include "IonNeutralAdvectionDiffusionTerm.tcc"
}

#endif/*_DREAM_EQUATION_ION_NEUTRAL_ADVECTION_DIFFUSION_TERM_HPP*/
