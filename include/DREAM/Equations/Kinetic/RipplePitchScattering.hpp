#ifndef _DREAM_EQUATIONS_KINETIC_RIPPLE_PITCH_SCATTERING_HPP
#define _DREAM_EQUATIONS_KINETIC_RIPPLE_PITCH_SCATTERING_HPP

#include "DREAM/MultiInterpolator1D.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class RipplePitchScattering : public FVM::DiffusionTerm {
    private:
        OptionConstants::eqterm_ripple_mode mode;
        real_t deltaCoils;              // Distance between toroidal field coils

        len_t nModes;
        const int_t *m, *n;             // Number of modes; poloidal mode numbers; toroidal mode numbers
        DREAM::MultiInterpolator1D *dB_B;    

        real_t **p_mn=nullptr;          // Resonant momentum (size nModes-by-nr)

        // If 'true', includes the poloidal field component when
        // evaluating the resonant momentum
        const bool INCLUDE_POLOIDAL_FIELD=false;
    public:
        RipplePitchScattering(
            FVM::Grid*, enum OptionConstants::eqterm_ripple_mode, 
            enum OptionConstants::momentumgrid_type,
            const real_t, const len_t, const int_t*, const int_t*,
            DREAM::MultiInterpolator1D*
        );
        RipplePitchScattering(
            FVM::Grid*, enum OptionConstants::eqterm_ripple_mode,
            enum OptionConstants::momentumgrid_type,
            const len_t, const len_t, const int_t*, const int_t*,
            DREAM::MultiInterpolator1D*
        );
        ~RipplePitchScattering();

        void Allocate();
        void CalculateResonantMomentum();

        const len_t GetNumberOfModes() { return this->nModes; }
        real_t **GetResonantMomentum() { return this->p_mn; }
        const int_t *GetPoloidalModeNumbers() { return this->m; }
        const int_t *GetToroidalModeNumbers() { return this->n; }

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_RIPPLE_PITCH_SCATTERING_HPP*/
