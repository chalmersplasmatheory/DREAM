#ifndef _DREAM_EQUATIONS_KINETIC_WAVE_PITCH_SCATTERING_HPP
#define _DREAM_EQUATIONS_KINETIC_WAVE_PITCH_SCATTERING_HPP

#include "DREAM/MultiInterpolator1D.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class WavePitchScattering : public FVM::DiffusionTerm {
    private:
        DREAM::MultiInterpolator1D *ppar_res; // resonant parallel momentum
        DREAM::MultiInterpolator1D *Delta_ppar_res; // resonance width in parallel momentum
        DREAM::MultiInterpolator1D *Dxx_int; // integrated xi-xi diffusion term
        const real_t t_start; // start time wave
        const real_t t_end; // end time wave
    public:
        WavePitchScattering(
            FVM::Grid*, enum OptionConstants::eqterm_wave_mode,
            enum OptionConstants::momentumgrid_type,
            DREAM::MultiInterpolator1D*, DREAM::MultiInterpolator1D*, DREAM::MultiInterpolator1D*
            const real_t, const real_t
        );
        WavePitchScattering(
            FVM::Grid*, enum OptionConstants::eqterm_ripple_mode,
            enum OptionConstants::momentumgrid_type,
            DREAM::MultiInterpolator1D*, DREAM::MultiInterpolator1D*, DREAM::MultiInterpolator1D*
            const real_t, const real_t
        );
        ~WavePitchScattering();

        void Allocate();
            
        // additional functions !!!
        const len_t GetResonantMomentum { return this->ppar_res;}

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_WAVE_PITCH_SCATTERING_HPP*/
