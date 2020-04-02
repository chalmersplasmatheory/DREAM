#ifndef _DREAM_EQUATIONS_KINETIC_PITCH_SCATTER_TERM_HPP
#define _DREAM_EQUATIONS_KINETIC_PITCH_SCATTER_TERM_HPP



#include "FVM/config.h"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/CollisionFrequencyCreator.hpp"


namespace DREAM {
    class PitchScatterTerm
        : public FVM::DiffusionTerm {
    private:
        enum SimulationGenerator::momentumgrid_type gridtype;
        CollisionFrequencyCreator *collFreqs;
        EquationSystem *eqSys;
        FVM::Grid *grid;
    public:
        PitchScatterTerm(FVM::Grid*,CollisionFrequencyCreator*,EquationSystem*,enum SimulationGenerator::momentumgrid_type);
        
        
        virtual void Rebuild();
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_PITCH_SCATTER_TERM_HPP*/


