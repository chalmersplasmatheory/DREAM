/**
 * Implementation of ion equations.
 */

#include "DREAM/Settings/SimulationGenerator.hpp"


using namespace DREAM;
using namespace std;


#define MODULENAME "ions"

/**
 * Define options for the ions.
 *
 * set: Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_Ions(Settings *s) {
    const len_t dims[1] = {0};

    s->DefineSetting(MODULENAME "/Z", "List of atomic charge numbers", 1, dims, (int_t*)nullptr);
    s->DefineSetting(MODULENAME "/types", "Method to use for determining ion charge distributions", 1, dims, (int_t*)nullptr);

    DefineDataIonR(MODULENAME, s, "density");
    DefineDataIonRT(MODULENAME, s, "prescribed");
}

/**
 * Construct the equation governing the evolution of the
 * ion densities.
 */
void SimulationGenerator::ConstructEquation_Ions(EquationSystem *eqsys, Settings *s) {
    Grid *fluidGrid = eqsys->GetFluidGrid();

    len_t nZ, ntypes;
    int_t *Z  = s->GetIntegerArray(MODULENAME "/Z", 1, &nZ);
    int_t *itypes = s->GetIntegerArray(MODULENAME "/types", 1, &ntypes);

    enum OptionConstants::ion_data_type *types = new enum OptionConstants::ion_data_type[ntypes];
    for (len_t i = 0; i < ntypes; i++)
        types[i] = (enum OptionConstants::ion_data_type)itypes[i];

    IonHandler *ih = new IonHandler(fluidGrid->GetRadialGrid(), unknowns, Z, nZ);

    // TODO Load ion data; initialize IonHandler
}

