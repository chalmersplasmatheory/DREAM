#ifndef _DREAM_EQUATIONS_KINETIC_ELECTRIC_FIELD_TERM_HPP
#define _DREAM_EQUATIONS_KINETIC_ELECTRIC_FIELD_TERM_HPP


#include "FVM/config.h"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"


namespace DREAM {
    class ElectricFieldTerm
        : public FVM::AdvectionTerm {
    private:
        enum SimulationGenerator::momentumgrid_type gridtype;
        EquationSystem *eqSys;
        FVM::Grid *grid;
        len_t id_Eterm;
    public:
        ElectricFieldTerm(FVM::Grid*,EquationSystem*,enum SimulationGenerator::momentumgrid_type);
        
        
        virtual void Rebuild(const real_t) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_ELECTRIC_FIELD_TERM_HPP*/


