#ifndef _DREAM_EQUATIONS_PARTICLE_SOURCE_TERM_HPP
#define _DREAM_EQUATIONS_PARTICLE_SOURCE_TERM_HPP

#include "DREAM/Equations/Kinetic/FluidKineticSourceTerm.hpp"
#include <limits>
namespace DREAM {
    class AvalancheSourceRP
        : public FluidKineticSourceTerm {
    private:
        len_t id_ntot;
        real_t preFactor;
        real_t pCutoff;
    protected:
        virtual real_t GetSourceFunction(len_t ir, len_t i, len_t j) override;
        virtual real_t GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId) override;
    public:
        AvalancheSourceRP(FVM::Grid*, FVM::UnknownQuantityHandler*, real_t);

        real_t EvaluateRPSource(len_t ir, len_t i, len_t j);
        real_t EvaluateTotalKnockOnNumber(len_t ir, real_t pLower, real_t pUpper=std::numeric_limits<real_t>::infinity());
    };
}

#endif/*_DREAM_EQUATIONS_PARTICLE_SOURCE_TERM_HPP*/


