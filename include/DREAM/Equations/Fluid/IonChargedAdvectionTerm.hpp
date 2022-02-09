#ifndef _DREAM_EQUATION_ION_CHARGED_ADVECTION_TERM_HPP
#define _DREAM_EQUATION_ION_CHARGED_ADVECTION_TERM_HPP

#include "DREAM/Equations/Fluid/IonChargedAdvectionDiffusionTerm.hpp"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/AdvectionInterpolationCoefficient.hpp"

namespace DREAM {
	class IonChargedAdvectionTerm : public IonChargedAdvectionDiffusionTerm<FVM::AdvectionTerm> {
	private:
    std::vector<FVM::AdvectionInterpolationCoefficient*> deltars;

	protected:
	
	    void AllocateInterpolationCoefficients();
	    void DeallocateInterpolationCoefficients(); 
	    void RebuildInterpolationCoefficients(FVM::UnknownQuantityHandler*, real_t**, real_t**, real_t**);
        void SetAdvectionInterpolationMethod(
                FVM::AdvectionInterpolationCoefficient::adv_interpolation intp,
                OptionConstants::adv_jacobian_mode jac_mode, 
                len_t id, real_t damping_factor=1.0 
            );
		
	public:
		IonChargedAdvectionTerm(FVM::Grid *g, IonHandler *ihdl, const len_t iIon, bool allocCoefficients,
		    FVM::AdvectionInterpolationCoefficient::adv_interpolation intp, OptionConstants::adv_jacobian_mode jac_mode, 
            len_t id, real_t damping_factor=1.0 );
		~IonChargedAdvectionTerm();
		
		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
		
		virtual bool SetCSJacobianBlock(
		    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *nions,
		    const len_t iIon, const len_t Z0, const len_t rOffset
		) override;
		virtual void SetCSMatrixElements(
		    FVM::Matrix *mat, real_t *rhs, const len_t iIon, const len_t Z0, const len_t rOffset
		) override;
		virtual void SetCSVectorElements(
		    real_t *vec, const real_t *nions, const len_t iIon, const len_t Z0, const len_t rOffset
		) override;
	};
}
#endif/*_DREAM_EQUATION_ION_CHARGED_ADVECTION_TERM_HPP*/
