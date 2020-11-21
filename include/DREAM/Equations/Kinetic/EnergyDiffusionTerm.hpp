#ifndef _DREAM_EQUATIONS_KINETIC_ENERGY_DIFFUSION_TERM_HPP
#define _DREAM_EQUATIONS_KINETIC_ENERGY_DIFFUSION_TERM_HPP

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
//#include "DREAM/Equations/ParallelDiffusionFrequency.hpp"
#include "FVM/config.h"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class EnergyDiffusionTerm
        : public FVM::DiffusionTerm {
    private:
        enum OptionConstants::momentumgrid_type gridtype;
        ParallelDiffusionFrequency *nuPar;
        virtual void SetPartialDiffusionTerm(len_t derivId, len_t nMultiples) override;
    public:
        EnergyDiffusionTerm(FVM::Grid*,CollisionQuantityHandler*,
            enum OptionConstants::momentumgrid_type, FVM::UnknownQuantityHandler*,bool);
        
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_ENERGY_DIFFUSION_TERM_HPP*/


