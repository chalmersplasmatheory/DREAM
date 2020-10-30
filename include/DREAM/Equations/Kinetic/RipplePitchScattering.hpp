#ifndef _DREAM_EQUATIONS_KINETIC_RIPPLE_PITCH_SCATTERING_HPP
#define _DREAM_EQUATIONS_KINETIC_RIPPLE_PITCH_SCATTERING_HPP

#include "DREAM/IonInterpolator1D.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class RipplePitchScattering : public FVM::DiffusionTerm {
    private:
        real_t deltaCoils;              // Distance between toroidal field coils

        len_t nModes;
        const int_t *m, *n;             // Number of modes; poloidal mode numbers; toroidal mode numbers
        DREAM::IonInterpolator1D *dB_B;    

        real_t **p_mn=nullptr;          // Resonant momentum (size nModes-by-nr)

        // If 'true', includes the poloidal field component when
        // evaluating the resonant momentum
        const bool INCLUDE_POLOIDAL_FIELD=false;

    public:
        RipplePitchScattering(
            FVM::Grid*, enum OptionConstants::momentumgrid_type,
            const real_t, const len_t, const int_t*, const int_t*,
            DREAM::IonInterpolator1D*
        );
        RipplePitchScattering(
            FVM::Grid*, enum OptionConstants::momentumgrid_type,
            const len_t, const len_t, const int_t*, const int_t*,
            DREAM::IonInterpolator1D*
        );
        ~RipplePitchScattering();

        void Allocate();
        void CalculateResonantMomentum();

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_RIPPLE_PITCH_SCATTERING_HPP*/
