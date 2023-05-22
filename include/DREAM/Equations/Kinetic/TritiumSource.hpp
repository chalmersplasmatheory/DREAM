#ifndef _DREAM_EQUATIONS_TRITIUM_SOURCE_HPP
#define _DREAM_EQUATIONS_TRITIUM_SOURCE_HPP

#include "DREAM/Equations/FluidSourceTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include <limits>

namespace DREAM {
    class TritiumSource
        : public FluidSourceTerm {
    
    public:
        enum SourceMode{
            SOURCE_MODE_FLUID,
            SOURCE_MODE_KINETIC
        };
    private:
        len_t id_nT;
        IonHandler *ions;
        real_t *sourceVec = nullptr;
        
        real_t *source;
        real_t pc, scaleFactor;
        
        SourceMode sourceMode;
	
    protected:
        static real_t integrand(real_t, void * );
        static real_t fbeta(real_t, void * );
        
        virtual real_t GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId) override;
    public:
        TritiumSource(FVM::Grid*, FVM::UnknownQuantityHandler*, real_t, real_t, SourceMode sm = SOURCE_MODE_KINETIC);
        
        virtual real_t GetSourceFunction(len_t ir, len_t i, len_t j) override;
        
        real_t EvaluateSource(len_t ir, len_t i, len_t j);
        static real_t EvaluateTotalTritiumNumber(real_t pLower, real_t pUpper=std::numeric_limits<real_t>::infinity());
        
        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual bool SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t* x) override;
    };
}

#endif/*_DREAM_EQUATIONS_TRITIUM_SOURCE_HPP*/
