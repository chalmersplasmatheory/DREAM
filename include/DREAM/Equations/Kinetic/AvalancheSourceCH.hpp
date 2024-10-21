#ifndef _DREAM_EQUATIONS_AVALANCHE_SOURCE_CH_HPP
#define _DREAM_EQUATIONS_AVALANCHE_SOURCE_CH_HPP

#include "DREAM/Equations/FluidSourceTerm.hpp"
#include <limits>
namespace DREAM {
    class AvalancheSourceCH
        : public FluidSourceTerm {
    public:
        enum CHSourceMode{
            CH_SOURCE_MODE_FLUID,
            CH_SOURCE_MODE_KINETIC
        };

        // Setting to determine whether to apply the source term for
        // negative or positive xi. In the "adaptive" case, this is determined
        // by the sign of the electric field.
        enum CHSourcePitchMode {
            CH_SOURCE_PITCH_ADAPTIVE,
            CH_SOURCE_PITCH_POSITIVE,
            CH_SOURCE_PITCH_NEGATIVE
        };
    private:
        len_t id_ntot;
        len_t id_Efield;
        real_t scaleFactor;
        real_t preFactor;
        real_t pCutoff;
        CHSourcePitchMode sourceXiMode;

    protected:
        virtual real_t GetSourceFunction(len_t ir, len_t i, len_t j) override;
        virtual real_t GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId) override;
    public:
        AvalancheSourceCH(FVM::Grid*, FVM::UnknownQuantityHandler*, real_t, real_t, CHSourceMode sm = CH_SOURCE_MODE_KINETIC, CHSourcePitchMode sem = CH_SOURCE_PITCH_ADAPTIVE);

        real_t EvaluateCHSource(len_t ir, len_t i, len_t j);
        static real_t EvaluateNormalizedTotalKnockOnNumber(real_t pLower, real_t pUpper=std::numeric_limits<real_t>::infinity());
    };
}

#endif/*_DREAM_EQUATIONS_AVALANCHE_SOURCE_CH_HPP*/


