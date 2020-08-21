/**
 * Implementation of various 'Simulation' helper routines.
 */

#include <string>
#include <softlib/SFile.h>
#include "DREAM/Simulation.hpp"


using namespace DREAM;
using namespace std;

/**
 * Constructor.
 */
Simulation::Simulation() {}


/**
 * Destructor.
 */
Simulation::~Simulation() {
    delete this->adas;
}

/**
 * Run this simulation.
 */
void Simulation::Run() {
    eqsys->Solve();
}

/**
 * Save the current state of this simulation to the
 * file with the previously specified name.
 */
void Simulation::Save() {
    this->Save(this->output_filename);
}

/**
 * Save the current state of this simulation.
 *
 * filename: Name of file to save simulation data to.
 */
void Simulation::Save(const string& filename) {
    SFile *sf = SFile::Create(filename, SFILE_MODE_WRITE);
    Save(sf);

    sf->Close();
}

/**
 * Save the current state of this simulation.
 *
 * sf: SFile object to use for saving the simulation.
 */
void Simulation::Save(SFile *sf) {
    // TODO Save settings

    // Save grids
    sf->CreateStruct("grid");
    eqsys->SaveGrids(sf, "/grid");
    
    // Save unknowns
    sf->CreateStruct("eqsys");
    eqsys->GetUnknownHandler()->SaveSFile(sf, "/eqsys", false);

    // Save ion metadata
    sf->CreateStruct("ionmeta");
    eqsys->SaveIonMetaData(sf, "/ionmeta");

    // Save timing information
    eqsys->SaveTimings(sf, "/timings");

    // Save "other" quantities (if any)
    OtherQuantityHandler *oqh = eqsys->GetOtherQuantityHandler();
    if (oqh->GetNRegistered() > 0) {
        sf->CreateStruct("other");
        oqh->SaveSFile(sf, "/other");
    }
}

