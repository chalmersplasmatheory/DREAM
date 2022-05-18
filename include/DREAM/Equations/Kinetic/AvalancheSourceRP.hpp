#ifndef _DREAM_EQUATIONS_AVALANCHE_SOURCE_RP_HPP
#define _DREAM_EQUATIONS_AVALANCHE_SOURCE_RP_HPP

#include "DREAM/Equations/FluidSourceTerm.hpp"
#include <limits>
namespace DREAM {
    class AvalancheSourceRP
        : public FluidSourceTerm {
    public:
        enum RPSourceMode{
            RP_SOURCE_MODE_FLUID,
            RP_SOURCE_MODE_KINETIC
        };

        // Setting to determine whether to apply the source term for
        // negative or positive xi. In the "adaptive" case, this is determined
        // by the sign of the electric field.
        enum RPSourcePitchMode {
            RP_SOURCE_PITCH_ADAPTIVE,
            RP_SOURCE_PITCH_POSITIVE,
            RP_SOURCE_PITCH_NEGATIVE
        };
    private:
        len_t id_ntot;
        len_t id_Efield;
        real_t scaleFactor;
        real_t preFactor;
        real_t pCutoff;
        RPSourceMode sourceMode;
        RPSourcePitchMode sourceXiMode;

    protected:
        virtual real_t GetSourceFunction(len_t ir, len_t i, len_t j) override;
        virtual real_t GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId) override;
    public:
        AvalancheSourceRP(FVM::Grid*, FVM::UnknownQuantityHandler*, real_t, real_t, RPSourceMode sm = RP_SOURCE_MODE_KINETIC, RPSourcePitchMode sem = RP_SOURCE_PITCH_ADAPTIVE);

        real_t EvaluateRPSource(len_t ir, len_t i, len_t j);
        static real_t EvaluateNormalizedTotalKnockOnNumber(real_t pLower, real_t pUpper=std::numeric_limits<real_t>::infinity());
    };
}

#endif/*_DREAM_EQUATIONS_AVALANCHE_SOURCE_RP_HPP*/


