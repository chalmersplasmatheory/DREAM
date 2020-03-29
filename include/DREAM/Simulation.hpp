#ifndef _DREAM_SIMULATION_HPP
#define _DREAM_SIMULATION_HPP
/**
 * This file contains the definition of the DREAM Simulation
 * class, which is the most central object in the code.
 * The simulation contains
 *
 *   - information about prescribed plasma parameters
 *   - the equation system to solve
 *   - the modules to use for advancing the system in time
 *   - ...
 */

#include "DREAM/EquationSystem.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAM {
    class Simulation {
    private:
        EquationSystem *eqsys;

    public:
        Simulation();
        ~Simulation();

        void Run();

        void SetEquationSystem(EquationSystem *e) { this->eqsys = eqsys; }
    };
}

#endif/*_DREAM_SIMULATION_HPP*/
