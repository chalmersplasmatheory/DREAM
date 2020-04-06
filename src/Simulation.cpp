/**
 * Implementation of various 'Simulation' helper routines.
 */

#include "DREAM/Simulation.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
Simulation::Simulation() {}

/**
 * Run this simulation.
 */
void Simulation::Run() {
    eqsys->ProcessSystem();
    eqsys->Solve();
}

