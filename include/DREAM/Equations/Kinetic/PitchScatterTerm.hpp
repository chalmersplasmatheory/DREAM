#ifndef _DREAM_EQUATIONS_KINETIC_PITCH_SCATTER_TERM_HPP
#define _DREAM_EQUATIONS_KINETIC_PITCH_SCATTER_TERM_HPP



#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "FVM/config.h"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
    class PitchScatterTerm
        : public FVM::DiffusionTerm {
    private:
        enum SimulationGenerator::momentumgrid_type gridtype;
        CollisionQuantityHandler *collQty;
        EquationSystem *eqSys;
        FVM::Grid *grid;
    public:
        PitchScatterTerm(FVM::Grid*,CollisionQuantityHandler*,EquationSystem*,enum SimulationGenerator::momentumgrid_type);
        
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_PITCH_SCATTER_TERM_HPP*/


