#ifndef _DREAM_EQUATIONS_KINETIC_RIPPLE_PITCH_SCATTERING_HPP
#define _DREAM_EQUATIONS_KINETIC_RIPPLE_PITCH_SCATTERING_HPP

#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class RipplePitchScattering : public FVM::DiffusionTerm {
    private:
        len_t nCoils;           // Number of toroidal field coils

        len_t nModes;
        const len_t *m, *n;             // Number of modes; poloidal mode numbers; toroidal mode numbers
        FVM::Interpolator1D **dB_B;     // Size

        real_t **p_mn=nullptr;          // Resonant momentum (size nModes-by-nr)

        // If 'true', includes the poloidal field component when
        // evaluating the resonant momentum
        const bool INCLUDE_POLOIDAL_FIELD=false;

    public:
        RipplePitchScattering(
            FVM::Grid*, const len_t, const len_t, const len_t*, const len_t*,
            FVM::Interpolator1D**
        );
        ~RipplePitchScattering();

        void Allocate();
        void CalculateResonantMomentum();

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_RIPPLE_PITCH_SCATTERING_HPP*/
