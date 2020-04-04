#ifndef _DREAM_EQUATIONS_KINETIC_SLOWING_DOWN_TERM_HPP
#define _DREAM_EQUATIONS_KINETIC_SLOWING_DOWN_TERM_HPP


#include "FVM/config.h"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"


namespace DREAM {
    class SlowingDownTerm
        : public FVM::AdvectionTerm {
    private:
        enum SimulationGenerator::momentumgrid_type gridtype;
        CollisionQuantityHandler *collQty;
    public:
        SlowingDownTerm(FVM::Grid*,CollisionQuantityHandler*,enum SimulationGenerator::momentumgrid_type);
        
        
        virtual void Rebuild();
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_SLOWING_DOWN_TERM_HPP*/


