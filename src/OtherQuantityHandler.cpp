/**
 * The 'OtherQuantityHandler' keeps track of the time evolution of various
 * quantities throughout the simulation. "Other" quantities are those which are
 * not solved for in the equation system and include things such as
 *
 *   - bounce coefficients
 *   - collision frequencies
 *   - operaetor coefficients
 *   - ...
 *
 * HOW TO ADD SUPPORT FOR A NEW QUANTITY:
 *   1. Give the quantity an ID by adding an entry to the
 *      'quantity_id' enum in the 'OtherQuantityHandler.hpp'
 *      file associated with this class.
 *   2. Add appropriate handlers to the three methods
 *        - _ConstructQuantity()
 *        - _GetName()
 *        - _StoreQuantity()
 *      which are declared further down in this file.
 */

#include <map>
#include "DREAM/OtherQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;
using FVM::QuantityData;
using namespace std;


/**
 * Constructor.
 */
OtherQuantityHandler::OtherQuantityHandler(
    CollisionQuantityHandler *cqtyHottail, CollisionQuantityHandler *cqtyRunaway,
    FVM::Grid *fluidGrid, FVM::Grid *hottailGrid, FVM::Grid *runawayGrid
) : cqtyHottail(cqtyHottail), cqtyRunaway(cqtyRunaway),
    fluidGrid(fluidGrid), hottailGrid(hottailGrid), runawayGrid(runawayGrid) {
}


/**
 * Destructor.
 */
OtherQuantityHandler::~OtherQuantityHandler() {
    for (auto it = this->data.begin(); it != this->data.end(); it++)
        delete it->second;
}

/**
 * Register a new "other" quantity to keep track of.
 *
 * id: Internal ID of quantity.
 */
void OtherQuantityHandler::RegisterQuantity(enum quantity_id id) {
    this->data[id] = _ConstructQuantity(id);
}

/**
 * Register all available "other" quantities.
 */
void OtherQuantityHandler::RegisterAllQuantities() {
    for (int i = 1; i != OTHER_QTY_LAST; i++) {
        enum quantity_id id = static_cast<enum quantity_id>(i);
        this->RegisterQuantity(id);
    }
}

/**
 * Store the values of all registered quantities in the
 * current time step.
 */
void OtherQuantityHandler::StoreAll(const real_t t) {
    for (auto it = this->data.begin(); it != this->data.end(); it++)
        this->_StoreQuantity(t, it->first, it->second);
}

/**
 * Save stored data to file.
 */
void OtherQuantityHandler::SaveSFile(SFile *sf, const std::string& path) {
    string group = path;
    if (path.back() != '/')
        group += '/';

    bool runawayCreated = false;
    bool hottailCreated = false;

    // Loop over and save stored quantities
    for (auto it = data.begin(); it != data.end(); it++) {
        const string name = _GetName(it->first);

        // Should we create any new groups first?
        if (!runawayCreated && name.substr(0, 8) == "runaway/")
            sf->CreateStruct(group+"runaway");
        else if (!hottailCreated && name.substr(0, 8) == "hottail/")
            sf->CreateStruct(group+"hottail");

        it->second->SaveSFile(sf, name, path);
    }
}

/********************************
 * IMPLEMENTATION OF QUANTITIES *
 ********************************/
/**
 * Construct a new 'QuantityData' object for the
 * specified quantity.
 */
QuantityData *OtherQuantityHandler::_ConstructQuantity(enum quantity_id id) {
    switch (id) {
        case OTHER_QTY_NU_S_HOTTAIL_FR:
            return new QuantityData(this->hottailGrid, 1, FVM::FLUXGRIDTYPE_RADIAL);
        case OTHER_QTY_NU_S_HOTTAIL_F1:
            return new QuantityData(this->hottailGrid, 1, FVM::FLUXGRIDTYPE_P1);
        case OTHER_QTY_NU_S_HOTTAIL_F2:
            return new QuantityData(this->hottailGrid, 1, FVM::FLUXGRIDTYPE_P2);
        case OTHER_QTY_NU_S_HOTTAIL:
            return new QuantityData(this->hottailGrid);

        case OTHER_QTY_NU_S_RUNAWAY_FR:
            return new QuantityData(this->runawayGrid, 1, FVM::FLUXGRIDTYPE_RADIAL);
        case OTHER_QTY_NU_S_RUNAWAY_F1:
            return new QuantityData(this->runawayGrid, 1, FVM::FLUXGRIDTYPE_P1);
        case OTHER_QTY_NU_S_RUNAWAY_F2:
            return new QuantityData(this->runawayGrid, 1, FVM::FLUXGRIDTYPE_P2);
        case OTHER_QTY_NU_S_RUNAWAY:
            return new QuantityData(this->runawayGrid);

        default:
            throw OtherQuantityException("Unrecognized other quantity ID: %d", id);
    }
}

/**
 * Returns the name of the specified "other" quantity.
 */
const char *OtherQuantityHandler::_GetName(enum quantity_id id) {
    switch (id) {
        case OTHER_QTY_NU_S_HOTTAIL_FR: return "hottail/nu_S_fr";
        case OTHER_QTY_NU_S_HOTTAIL_F1: return "hottail/nu_S_f1";
        case OTHER_QTY_NU_S_HOTTAIL_F2: return "hottail/nu_S_f2";
        case OTHER_QTY_NU_S_HOTTAIL:    return "hottail/nu_S";

        case OTHER_QTY_NU_S_RUNAWAY_FR: return "runaway/nu_S_fr";
        case OTHER_QTY_NU_S_RUNAWAY_F1: return "runaway/nu_S_f1";
        case OTHER_QTY_NU_S_RUNAWAY_F2: return "runaway/nu_S_f2";
        case OTHER_QTY_NU_S_RUNAWAY:    return "runaway/nu_S";

        default:
            throw OtherQuantityException("Unrecognized other quantity ID: %d", id);
    }
}

/**
 * Store data for the specified quantity in the current
 * time step.
 */
void OtherQuantityHandler::_StoreQuantity(
    const real_t t, enum quantity_id id, QuantityData *qd
) {
    // XXX here we assume that all momentum grids are the same
    const len_t nr_ht = this->hottailGrid->GetNr();
    const len_t n1_ht = this->hottailGrid->GetMomentumGrid(0)->GetNp1();
    const len_t n2_ht = this->hottailGrid->GetMomentumGrid(0)->GetNp2();

    const len_t nr_re = this->runawayGrid->GetNr();
    const len_t n1_re = this->runawayGrid->GetMomentumGrid(0)->GetNp1();
    const len_t n2_re = this->runawayGrid->GetMomentumGrid(0)->GetNp2();

    switch (id) {
        // Slowing-down collision frequency on hot-tail grid
        case OTHER_QTY_NU_S_HOTTAIL_FR:
            qd->Store(nr_ht+1, n1_ht*n2_ht, this->cqtyHottail->GetNuS()->GetValue_fr());
            break;
        case OTHER_QTY_NU_S_HOTTAIL_F1:
            qd->Store(nr_ht, (n1_ht+1)*n2_ht, this->cqtyHottail->GetNuS()->GetValue_f1());
            break;
        case OTHER_QTY_NU_S_HOTTAIL_F2:
            qd->Store(nr_ht, n1_ht*(n2_ht+1), this->cqtyHottail->GetNuS()->GetValue_f2());
            break;
        case OTHER_QTY_NU_S_HOTTAIL:
            qd->Store(nr_ht, n1_ht*n2_ht, this->cqtyHottail->GetNuS()->GetValue());
            break;

        // Slowing-down collision frequency on runaway grid
        case OTHER_QTY_NU_S_RUNAWAY_FR:
            qd->Store(nr_re+1, n1_re*n2_re, this->cqtyRunaway->GetNuS()->GetValue_fr());
            break;
        case OTHER_QTY_NU_S_RUNAWAY_F1:
            qd->Store(nr_re, (n1_re+1)*n2_re, this->cqtyRunaway->GetNuS()->GetValue_f1());
            break;
        case OTHER_QTY_NU_S_RUNAWAY_F2:
            qd->Store(nr_re, n1_re*(n2_re+1), this->cqtyRunaway->GetNuS()->GetValue_f2());
            break;
        case OTHER_QTY_NU_S_RUNAWAY:
            qd->Store(nr_re, n1_re*n2_re, this->cqtyRunaway->GetNuS()->GetValue());
            break;

        default:
            throw OtherQuantityException("Unrecognized other quantity ID: %d", id);
    }

    qd->SaveStep(t);
}
