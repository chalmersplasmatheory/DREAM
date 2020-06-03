/**
 * The 'OtherQuantityHandler' keeps track of the time evolution of various
 * quantities throughout the simulation. "Other" quantities are those which are
 * not solved for in the equation system and include things such as
 *
 *   - bounce coefficients
 *   - collision frequencies
 *   - operator coefficients
 *   - ...
 *
 * HOW TO ADD SUPPORT FOR A NEW QUANTITY:
 *   1. Add a handler for the quantity to the "DefineQuantities()" method at
 *      the end of this file.
 *   2. Add the name of the quantity to the 'QUANTITIES' list in
 *      'py/DREAM/Settings/OtherQuantities.py' in the Python interface. This
 *      causes the interface to emit a warning if the user tries to add a
 *      non-existant quantity to the include.
 *
 * If the variable contains sections (i.e. the name contains forward slashes),
 * then you might also have to add code to 'SaveSFile()' for creating the
 * appropriate structs in the output file.
 */

#include <map>
#include "DREAM/OtherQuantity.hpp"
#include "DREAM/OtherQuantityHandler.hpp"
#include "DREAM/PostProcessor.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;
using FVM::QuantityData;
using namespace std;


/**
 * Constructor.
 */
OtherQuantityHandler::OtherQuantityHandler(
    CollisionQuantityHandler *cqtyHottail, CollisionQuantityHandler *cqtyRunaway,
    PostProcessor *postProcessor, RunawayFluid *REFluid,
    FVM::Grid *fluidGrid, FVM::Grid *hottailGrid, FVM::Grid *runawayGrid
) : cqtyHottail(cqtyHottail), cqtyRunaway(cqtyRunaway),
    postProcessor(postProcessor), REFluid(REFluid),
    fluidGrid(fluidGrid), hottailGrid(hottailGrid), runawayGrid(runawayGrid) {

    this->DefineQuantities();
}


/**
 * Destructor.
 */
OtherQuantityHandler::~OtherQuantityHandler() {
    for (auto it = this->all_quantities.begin(); it != this->all_quantities.end(); it++)
        delete *it;
}

/**
 * Returns an OtherQuantity by name.
 *
 * name: Name of quantity to return.
 */
OtherQuantity *OtherQuantityHandler::GetByName(const std::string& name) {
    for (auto it = this->all_quantities.begin(); it != this->all_quantities.end(); it++) {
        if ((*it)->GetName() == name)
            return *it;
    }

    return nullptr;
}

/**
 * Register a group of parameters.
 *
 * name: Name of group to register parameters from.
 */
bool OtherQuantityHandler::RegisterGroup(const std::string& name) {
    if (groups.find(name) != groups.end()) {
        vector<string>& grp = groups[name];
        for (auto it = grp.begin(); it != grp.end(); it++)
            this->RegisterQuantity(*it);
        return true;
    } else return false;
}

/**
 * Register a new "other" quantity to keep track of.
 *
 * name: Name of quantity to register.
 */
void OtherQuantityHandler::RegisterQuantity(const std::string& name) {
    OtherQuantity *oq = GetByName(name);

    if (oq == nullptr) {
        // Is the given name a group of parameters?
        // Try to register it!
        if (!RegisterGroup(name))
            throw OtherQuantityException("Unrecognized other quantity: '%s'.", name.c_str());
    // Skip the parameter if the grid is disabled
    } else if (oq->GetGrid() != nullptr)
        RegisterQuantity(oq);
}
void OtherQuantityHandler::RegisterQuantity(OtherQuantity *oq) {
    if (!oq->IsActive()) {
        oq->Activate();
        this->registered.push_back(oq);
    }
}

/**
 * Register all available "other" quantities.
 */
void OtherQuantityHandler::RegisterAllQuantities() {
    for (auto it = this->all_quantities.begin(); it != this->all_quantities.end(); it++)
        RegisterQuantity(*it);
}

/**
 * Store the values of all registered quantities in the
 * current time step.
 *
 * t: Current (simulation) time.
 */
void OtherQuantityHandler::StoreAll(const real_t t) {
    for (auto it = this->registered.begin(); it != this->registered.end(); it++)
        (*it)->Store(t);
}

/**
 * Save stored data to file.
 *
 * sf:   SFile object to save data to.
 * path: Path in SFile object to save data to (default: "").
 */
void OtherQuantityHandler::SaveSFile(SFile *sf, const std::string& path) {
    string group = path;
    if (path.back() != '/')
        group += '/';

    bool fluidCreated = false;
    bool hottailCreated = false;
    bool runawayCreated = false;

    // Loop over and save stored quantities
    for (auto it = this->registered.begin(); it != this->registered.end(); it++) {
        OtherQuantity *oq = *it;

        // Should we create any new groups first?
        if (!runawayCreated && oq->GetName().substr(0, 8) == "runaway/") {
            sf->CreateStruct(group+"runaway");
            runawayCreated = true;
        } else if (!hottailCreated && oq->GetName().substr(0, 8) == "hottail/") {
            sf->CreateStruct(group+"hottail");
            hottailCreated = true;
        } else if (!fluidCreated && oq->GetName().substr(0, 6) == "fluid/") {
            sf->CreateStruct(group+"fluid");
            fluidCreated = true;
        }
        oq->SaveSFile(sf, path);
    }
}

/********************************
 * IMPLEMENTATION OF QUANTITIES *
 ********************************/
/**
 * Define all other quantities.
 */
void OtherQuantityHandler::DefineQuantities() {
    // XXX here we assume that all momentum grids are the same
    const len_t nr_ht = (this->hottailGrid==nullptr ? 0 : this->hottailGrid->GetNr());
    const len_t n1_ht = (this->hottailGrid==nullptr ? 0 : this->hottailGrid->GetMomentumGrid(0)->GetNp1());
    const len_t n2_ht = (this->hottailGrid==nullptr ? 0 : this->hottailGrid->GetMomentumGrid(0)->GetNp2());

    const len_t nr_re = (this->runawayGrid==nullptr ? 0 : this->runawayGrid->GetNr());
    const len_t n1_re = (this->runawayGrid==nullptr ? 0 : this->runawayGrid->GetMomentumGrid(0)->GetNp1());
    const len_t n2_re = (this->runawayGrid==nullptr ? 0 : this->runawayGrid->GetMomentumGrid(0)->GetNp2());

    // HELPER MACROS (to make definitions more compact)
    // Define on fluid grid
    #define DEF_FL(NAME, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), fluidGrid, 1, FVM::FLUXGRIDTYPE_DISTRIBUTION, [this](QuantityData *qd) {FUNC}));
    #define DEF_FL_FR(NAME, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), fluidGrid, 1, FVM::FLUXGRIDTYPE_RADIAL, [this](QuantityData *qd) {FUNC}));

    // Define on hot-tail grid
    #define DEF_HT(NAME, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), hottailGrid, 1, FVM::FLUXGRIDTYPE_DISTRIBUTION, [this,nr_ht,n1_ht,n2_ht](QuantityData *qd) {FUNC}));
    #define DEF_HT_FR(NAME, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), hottailGrid, 1, FVM::FLUXGRIDTYPE_RADIAL, [this,nr_ht,n1_ht,n2_ht](QuantityData *qd) {FUNC}));
    #define DEF_HT_F1(NAME, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), hottailGrid, 1, FVM::FLUXGRIDTYPE_P1, [this,nr_ht,n1_ht,n2_ht](QuantityData *qd) {FUNC}));
    #define DEF_HT_F2(NAME, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), hottailGrid, 1, FVM::FLUXGRIDTYPE_P2, [this,nr_ht,n1_ht,n2_ht](QuantityData *qd) {FUNC}));

    // Define on runaway grid
    #define DEF_RE(NAME, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), runawayGrid, 1, FVM::FLUXGRIDTYPE_DISTRIBUTION, [this,nr_re,n1_re,n2_re](QuantityData *qd) {FUNC}));
    #define DEF_RE_FR(NAME, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), runawayGrid, 1, FVM::FLUXGRIDTYPE_RADIAL, [this,nr_re,n1_re,n2_re](QuantityData *qd) {FUNC}));
    #define DEF_RE_F1(NAME, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), runawayGrid, 1, FVM::FLUXGRIDTYPE_P1, [this,nr_re,n1_re,n2_re](QuantityData *qd) {FUNC}));
    #define DEF_RE_F2(NAME, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), runawayGrid, 1, FVM::FLUXGRIDTYPE_P2, [this,nr_re,n1_re,n2_re](QuantityData *qd) {FUNC}));
    
    // fluid/...
    DEF_FL("fluid/Eceff", qd->Store(this->REFluid->GetEffectiveCriticalField()););
    DEF_FL("fluid/EDreic", qd->Store(this->REFluid->GetDreicerElectricField()););
    DEF_FL("fluid/Ectot", qd->Store(this->REFluid->GetConnorHastieField_NOSCREENING()););
    DEF_FL("fluid/Ecfree", qd->Store(this->REFluid->GetConnorHastieField_COMPLETESCREENING()););
    DEF_FL("fluid/tauEERel", qd->Store(this->REFluid->GetElectronCollisionTimeRelativistic()););
    DEF_FL("fluid/tauEETh", qd->Store(this->REFluid->GetElectronCollisionTimeThermal()););
    DEF_FL("fluid/GammaAva", qd->Store(this->REFluid->GetAvalancheGrowthRate()););
    DEF_FL("fluid/lnLambdaC", qd->Store(this->REFluid->GetLnLambda()->GetLnLambdaC()););
    DEF_FL("fluid/lnLambdaT", qd->Store(this->REFluid->GetLnLambda()->GetLnLambdaT()););
    DEF_FL("fluid/runawayRate", qd->Store(this->postProcessor->GetRunawayRate()););

    // hottail/nu_s
    DEF_HT_F1("hottail/nu_s_f1", qd->Store(nr_ht,   (n1_ht+1)*n2_ht, this->cqtyHottail->GetNuS()->GetValue_f1()););
    DEF_HT_F2("hottail/nu_s_f2", qd->Store(nr_ht,   n1_ht*(n2_ht+1), this->cqtyHottail->GetNuS()->GetValue_f2()););
    DEF_HT_F1("hottail/nu_D_f1", qd->Store(nr_ht,   (n1_ht+1)*n2_ht, this->cqtyHottail->GetNuD()->GetValue_f1()););
    DEF_HT_F2("hottail/nu_D_f2", qd->Store(nr_ht,   n1_ht*(n2_ht+1), this->cqtyHottail->GetNuD()->GetValue_f2()););

    // runaway/nu_D
    DEF_RE_F1("runaway/nu_s_f1", qd->Store(nr_re,   (n1_re+1)*n2_re, this->cqtyRunaway->GetNuS()->GetValue_f1()););
    DEF_RE_F2("runaway/nu_s_f2", qd->Store(nr_re,   n1_re*(n2_re+1), this->cqtyRunaway->GetNuS()->GetValue_f2()););
    DEF_RE_F1("runaway/nu_D_f1", qd->Store(nr_re,   (n1_re+1)*n2_re, this->cqtyRunaway->GetNuD()->GetValue_f1()););
    DEF_RE_F2("runaway/nu_D_f2", qd->Store(nr_re,   n1_re*(n2_re+1), this->cqtyRunaway->GetNuD()->GetValue_f2()););


    // Declare groups of parameters (for registering
    // multiple parameters in one go)
    this->groups["fluid"] = {
        "fluid/Eceff", "fluid/GammaAva",
        "fluid/lnLambdaC", "fluid/lnLambdaT",
        "fluid/runawayRate", "fluid/Ectot",
        "fluid/Ecfree", "fluid/tauEERel",
        "fluid/tauEETh", "fluid/EDreic"
    };
    
    this->groups["nu_s"] = {
        "hottail/nu_s_f1", "hottail/nu_s_f2",
        "runaway/nu_s_f1", "runaway/nu_s_f2"
    };
    this->groups["nu_D"] = {
        "hottail/nu_D_f1", "hottail/nu_D_f2",
        "runaway/nu_D_f1", "runaway/nu_D_f2"
    };
}

