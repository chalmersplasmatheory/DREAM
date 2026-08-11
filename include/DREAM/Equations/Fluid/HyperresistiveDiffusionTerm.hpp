#ifndef _DREAM_EQUATIONS_POLOIDAL_FLUX_HYPERRESISTIVE_DIFFUSION_TERM_HPP
#define _DREAM_EQUATIONS_POLOIDAL_FLUX_HYPERRESISTIVE_DIFFUSION_TERM_HPP

#include "FVM/config.h"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"

namespace DREAM {
    class HyperresistiveDiffusionTerm
        : public FVM::DiffusionTerm {
    private:
		IonHandler *ions;
        FVM::Interpolator1D *Lambda=nullptr; 
		FVM::Interpolator1D *dBB=nullptr;
		
		real_t *lambda_buf = nullptr;
	
	protected:
		virtual const real_t *EvaluateLambda(const real_t t);

		void BuildCoefficient(const real_t*, real_t**);

    public:
        HyperresistiveDiffusionTerm(FVM::Grid*, IonHandler*, FVM::Interpolator1D*, FVM::Interpolator1D*);
		virtual ~HyperresistiveDiffusionTerm();
        
        virtual const real_t *GetLambda() const;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_POLOIDAL_FLUX_HYPERRESISTIVE_DIFFUSION_TERM_HPP*/


