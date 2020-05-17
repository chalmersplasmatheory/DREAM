#ifndef _DREAM_EQUATIONS_KINETIC_ELECTRIC_FIELD_DIFFUSION_TERM_HPP
#define _DREAM_EQUATIONS_KINETIC_ELECTRIC_FIELD_DIFFUSION_TERM_HPP


#include "FVM/config.h"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"


namespace DREAM {
    class ElectricFieldDiffusionTerm
        : public FVM::DiffusionTerm {
    private:
        PitchScatterFrequency *nuD;
        FVM::Grid *grid;
        len_t id_Eterm;

    public:
        ElectricFieldDiffusionTerm(FVM::Grid*, CollisionQuantityHandler*, FVM::UnknownQuantityHandler*);
        
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_ELECTRIC_FIELD_DIFFUSION_TERM_HPP*/


