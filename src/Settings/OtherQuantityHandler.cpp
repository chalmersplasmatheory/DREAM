/**
 * Construct an 'OtherQuantityHandler' object.
 */

#include "DREAM/OtherQuantityHandler.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Settings/Settings.hpp"


#define MODULENAME "other"

using namespace DREAM;

/**
 * Define all options available for the 'OtherQuantityHandler'.
 */
void SimulationGenerator::DefineOptions_OtherQuantities(Settings *s) {
    const len_t dims[1] = {0};

    s->DefineSetting(MODULENAME "/include", "List of IDs of other quantities to include", 1, dims, (int_t*)nullptr);
}

/**
 * Construct an 'OtherQuantityHandler' object.
 */
void SimulationGenerator::ConstructOtherQuantityHandler(
    EquationSystem *eqsys, Settings *s
) {
    OtherQuantityHandler *oqh = new OtherQuantityHandler(
        eqsys->GetHotTailCollisionHandler(), eqsys->GetRunawayCollisionHandler(),
        eqsys->GetFluidGrid(), eqsys->GetHotTailGrid(), eqsys->GetRunawayGrid()
    );

    len_t nother = 0;
    const int_t *other = s->GetIntegerArray(MODULENAME "/include", 1, &nother);

    // Catch-all (register all quantities)
    if (nother == 1 && *other < 0) {
        oqh->RegisterAllQuantities();
    } else if (nother > 0) {
        for (len_t i = 0; i < nother; i++)
            oqh->RegisterQuantity(static_cast<enum OtherQuantityHandler::quantity_id>(other[i]));
    }
    
    eqsys->SetOtherQuantityHandler(oqh);
}

