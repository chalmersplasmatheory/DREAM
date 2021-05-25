#ifndef _DREAM_EQUATION_ION_CHARGED_ADVECTION_DIFFUSION_TERM_HPP
#define _DREAM_EQUATION_ION_CHARGED_ADVECTION_DIFFUSION_TERM_HPP

#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"

namespace DREAM {
	template<class T>
    class IonChargedAdvectionDiffusionTerm : public IonEquationTerm<T> {
    protected:
		real_t **CoeffsAllCS = nullptr;
		len_t Z0ForPartials;
		
		void Allocate();
		void Deallocate();
		virtual void SetCoeffs(const len_t Z0)=0;
		virtual void SetCoeffsAllCS(const real_t t)=0;
		virtual void SetDiffCoeffsAllCS(const real_t t)=0;
		
		
	public:
		IonChargedAdvectionDiffusionTerm(FVM::Grid *g, IonHandler *ihdl, const len_t iIon);
		virtual ~IonChargedAdvectionDiffusionTerm();
		
		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
		
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
	}
	
	#include "IonChargedAdvectionDiffusionTerm.tcc"
}

#endif/*_DREAM_EQUATION_ION_CHARGED_ADVECTION_DIFFUSION_TERM_HPP*/
