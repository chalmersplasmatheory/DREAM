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

namespace DREAM { class Simulation; }

#include <string>
#include <softlib/SFile.h>
#include "DREAM/ADAS.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/NIST.hpp"
#include "DREAM/OutputGenerator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAM {
    class Simulation {
    private:
        ADAS *adas;
        NIST *nist;
        EquationSystem *eqsys;
		OutputGenerator *outgen=nullptr;

    public:
        Simulation();
        ~Simulation();

        void Run();

        ADAS *GetADAS() { return this->adas; }
        NIST *GetNIST() { return this->nist; }
        EquationSystem *GetEquationSystem() { return this->eqsys; }

        void SetADAS(ADAS *a) { this->adas = a; }
        void SetNIST(NIST *n) { this->nist = n; }
        void SetEquationSystem(EquationSystem *e) { this->eqsys = e; }

        void Save();
        void SetOutputGenerator(OutputGenerator*);
    };
}

#endif/*_DREAM_SIMULATION_HPP*/
