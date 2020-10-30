/**
 * Routines for saving things from the EquationSystem object.
 */

#include <string>
#include <softlib/SFile.h>
#include "DREAM/EquationSystem.hpp"


using namespace DREAM;
using namespace std;


/**
 * Save timinig information from the simulation.
 *
 * sf:   SFile object to saving timing information to.
 * name: Name of group to store information in.
 */
void EquationSystem::SaveTimings(SFile *sf, const string& name) {
    if (!this->timingFile) return;

    sf->CreateStruct(name);

    // Total simulation time
    sf->WriteScalar(name+"/total", this->simulationTime);
    sf->WriteAttribute_string(name+"/total", "desc", "Total simulation time");

    string path = name + "/solver";
    sf->CreateStruct(path);
    this->solver->SaveTimings(sf, path);

    path = name + "/runawayfluid";
    sf->CreateStruct(path);
    this->REFluid->SaveTimings(sf, path);
}

