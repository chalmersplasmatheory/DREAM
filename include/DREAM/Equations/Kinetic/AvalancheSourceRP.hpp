#ifndef _DREAM_EQUATIONS_PARTICLE_SOURCE_TERM_HPP
#define _DREAM_EQUATIONS_PARTICLE_SOURCE_TERM_HPP

#include "DREAM/Equations/Kinetic/FluidKineticSourceTerm.hpp"

namespace DREAM {
    class AvalancheSourceRP
        : public FluidKineticSourceTerm {
    private:
        len_t id_ntot;
        real_t EvaluateRPSource(len_t ir, len_t i, len_t j);
    protected:
        virtual real_t GetSourceFunction(len_t ir, len_t i, len_t j) override;
        virtual real_t GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId) override;
    public:
        AvalancheSourceRP(FVM::Grid*, FVM::UnknownQuantityHandler*);
        
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual void SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t* x) override;
    };
}

#endif/*_DREAM_EQUATIONS_PARTICLE_SOURCE_TERM_HPP*/


