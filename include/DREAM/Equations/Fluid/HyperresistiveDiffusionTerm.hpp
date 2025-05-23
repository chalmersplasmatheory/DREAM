#ifndef _DREAM_EQUATIONS_POLOIDAL_FLUX_HYPERRESISTIVE_DIFFUSION_TERM_HPP
#define _DREAM_EQUATIONS_POLOIDAL_FLUX_HYPERRESISTIVE_DIFFUSION_TERM_HPP

#include "FVM/config.h"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class HyperresistiveDiffusionTerm
        : public FVM::DiffusionTerm {
    private:
        FVM::Interpolator1D *Lambda=nullptr; 
	
	protected:
		virtual const real_t *EvaluateLambda(const real_t t) { return this->Lambda->Eval(t); }

		void BuildCoefficient(const real_t*, real_t**);

    public:
        HyperresistiveDiffusionTerm(FVM::Grid*, FVM::Interpolator1D*);
		virtual ~HyperresistiveDiffusionTerm();
        
        virtual const real_t *GetLambda() const { return Lambda->GetBuffer(); }
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_POLOIDAL_FLUX_HYPERRESISTIVE_DIFFUSION_TERM_HPP*/


