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
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/PostProcessor.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/OptionConstants.hpp"


using namespace DREAM;
using FVM::QuantityData;
using namespace std;


/**
 * Constructor.
 */
OtherQuantityHandler::OtherQuantityHandler(
    CollisionQuantityHandler *cqtyHottail, CollisionQuantityHandler *cqtyRunaway,
    PostProcessor *postProcessor, RunawayFluid *REFluid, FVM::UnknownQuantityHandler *unknowns, std::vector<UnknownQuantityEquation*> *unknown_equations, Settings *s,
    FVM::Grid *fluidGrid, FVM::Grid *hottailGrid, FVM::Grid *runawayGrid
) : cqtyHottail(cqtyHottail), cqtyRunaway(cqtyRunaway),
    postProcessor(postProcessor), REFluid(REFluid), unknowns(unknowns), unknown_equations(unknown_equations), s(s),
    fluidGrid(fluidGrid), hottailGrid(hottailGrid), runawayGrid(runawayGrid) {

    id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    eqn_Tcold = unknown_equations->at(id_Tcold);
    id_term_rad=0;

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

    const real_t *x = unknowns->GetUnknownData(id_ncold);

    // HELPER MACROS (to make definitions more compact)
    // Define on fluid grid
    #define DEF_FL(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), fluidGrid, 1, FVM::FLUXGRIDTYPE_DISTRIBUTION, [this,x](QuantityData *qd) {FUNC}));
    #define DEF_FL_FR(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), fluidGrid, 1, FVM::FLUXGRIDTYPE_RADIAL, [this](QuantityData *qd) {FUNC}));

    // Define on hot-tail grid
    #define DEF_HT(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), hottailGrid, 1, FVM::FLUXGRIDTYPE_DISTRIBUTION, [this,nr_ht,n1_ht,n2_ht](QuantityData *qd) {FUNC}));
    #define DEF_HT_FR(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), hottailGrid, 1, FVM::FLUXGRIDTYPE_RADIAL, [this,nr_ht,n1_ht,n2_ht](QuantityData *qd) {FUNC}));
    #define DEF_HT_F1(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), hottailGrid, 1, FVM::FLUXGRIDTYPE_P1, [this,nr_ht,n1_ht,n2_ht](QuantityData *qd) {FUNC}));
    #define DEF_HT_F2(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), hottailGrid, 1, FVM::FLUXGRIDTYPE_P2, [this,nr_ht,n1_ht,n2_ht](QuantityData *qd) {FUNC}));

    // Define on runaway grid
    #define DEF_RE(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), runawayGrid, 1, FVM::FLUXGRIDTYPE_DISTRIBUTION, [this,nr_re,n1_re,n2_re](QuantityData *qd) {FUNC}));
    #define DEF_RE_FR(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), runawayGrid, 1, FVM::FLUXGRIDTYPE_RADIAL, [this,nr_re,n1_re,n2_re](QuantityData *qd) {FUNC}));
    #define DEF_RE_F1(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), runawayGrid, 1, FVM::FLUXGRIDTYPE_P1, [this,nr_re,n1_re,n2_re](QuantityData *qd) {FUNC}));
    #define DEF_RE_F2(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), runawayGrid, 1, FVM::FLUXGRIDTYPE_P2, [this,nr_re,n1_re,n2_re](QuantityData *qd) {FUNC}));
    
    // fluid/...
    DEF_FL("fluid/Eceff", "Effective critical electric field [V/m]", qd->Store(this->REFluid->GetEffectiveCriticalField()););
    DEF_FL("fluid/Ecfree", "Connor-Hastie threshold field (calculated with n=n_free) [V/m]", qd->Store(this->REFluid->GetConnorHastieField_COMPLETESCREENING()););
    DEF_FL("fluid/Ectot", "Connor-Hastie threshold field (calculated with n=n_tot) [V/m]", qd->Store(this->REFluid->GetConnorHastieField_NOSCREENING()););
    DEF_FL("fluid/EDreic", "Dreicer electric field [V/m]", qd->Store(this->REFluid->GetDreicerElectricField()););
    DEF_FL("fluid/GammaAva", "Avalanche growth rate [s^-1]", qd->Store(this->REFluid->GetAvalancheGrowthRate()););
    DEF_FL("fluid/GammaDreicer", "Dreicer runaway rate [s^-1]", qd->Store(this->REFluid->GetDreicerRunawayRate()););
    DEF_FL("fluid/GammaCompton", "Compton runaway rate [s^-1]", qd->Store(this->REFluid->GetComptonRunawayRate()););
    DEF_FL("fluid/pCrit", "Critical momentum for avalanche (in units of mc)", qd->Store(this->REFluid->GetEffectiveCriticalRunawayMomentum()););
    DEF_FL("fluid/lnLambdaC", "Coulomb logarithm (relativistic)", qd->Store(this->REFluid->GetLnLambda()->GetLnLambdaC()););
    DEF_FL("fluid/lnLambdaT", "Coulomb logarithm (thermal)", qd->Store(this->REFluid->GetLnLambda()->GetLnLambdaT()););
    DEF_FL("fluid/runawayRate", "Total runaway rate, dn_RE / dt", qd->Store(this->postProcessor->GetRunawayRate()););
    DEF_FL("fluid/tauEERel", "Relativistic electron collision time (4*pi*lnL*n_cold*r^2*c)^-1 [s]", qd->Store(this->REFluid->GetElectronCollisionTimeRelativistic()););
    DEF_FL("fluid/tauEETh", "Thermal electron collision time (tauEERel * [2T/mc^2]^1.5) [s]", qd->Store(this->REFluid->GetElectronCollisionTimeThermal()););
    DEF_FL("fluid/conductivity", "Electric conductivity in SI, Sauter formula (based on Braams)", qd->Store(this->REFluid->GetElectricConductivity()););
    DEF_FL("fluid/Zeff", "Effective charge", qd->Store(this->REFluid->GetIonHandler()->evaluateZeff()););

    enum OptionConstants::uqty_T_cold_eqn type = (enum OptionConstants::uqty_T_cold_eqn)s->GetInteger("eqsys/T_cold/type");
    if (type==OptionConstants::UQTY_T_COLD_SELF_CONSISTENT)
        DEF_FL("fluid/radiation", "Radiated power density", qd->Store(this->eqn_Tcold->GetEquation(id_ncold)->GetVectorElementsSingleEquationTerm(id_term_rad,x)););

    // hottail/...
    DEF_HT_F1("hottail/nu_s_f1", "Slowing down frequency (on p1 flux grid) [s^-1]", qd->Store(nr_ht,   (n1_ht+1)*n2_ht, this->cqtyHottail->GetNuS()->GetValue_f1()););
    DEF_HT_F2("hottail/nu_s_f2", "Slowing down frequency (on p2 flux grid) [s^-1]", qd->Store(nr_ht,   n1_ht*(n2_ht+1), this->cqtyHottail->GetNuS()->GetValue_f2()););
    DEF_HT_F1("hottail/nu_D_f1", "Pitch-angle scattering frequency (on p1 flux grid) [s^-1]", qd->Store(nr_ht,   (n1_ht+1)*n2_ht, this->cqtyHottail->GetNuD()->GetValue_f1()););
    DEF_HT_F2("hottail/nu_D_f2", "Pitch-angle scattering frequency (on p2 flux grid) [s^-1]", qd->Store(nr_ht,   n1_ht*(n2_ht+1), this->cqtyHottail->GetNuD()->GetValue_f2()););
    DEF_HT_F1("hottail/lnLambda_ee_f1", "Coulomb logarithm for e-e collisions (on p1 flux grid)", qd->Store(nr_ht,   (n1_ht+1)*n2_ht, this->cqtyHottail->GetLnLambdaEE()->GetValue_f1()););
    DEF_HT_F2("hottail/lnLambda_ee_f2", "Coulomb logarithm for e-e collisions (on p2 flux grid)", qd->Store(nr_ht,   n1_ht*(n2_ht+1), this->cqtyHottail->GetLnLambdaEE()->GetValue_f2()););
    DEF_HT_F1("hottail/lnLambda_ei_f1", "Coulomb logarithm for e-i collisions (on p1 flux grid)", qd->Store(nr_ht,   (n1_ht+1)*n2_ht, this->cqtyHottail->GetLnLambdaEI()->GetValue_f1()););
    DEF_HT_F2("hottail/lnLambda_ei_f2", "Coulomb logarithm for e-i collisions (on p2 flux grid)", qd->Store(nr_ht,   n1_ht*(n2_ht+1), this->cqtyHottail->GetLnLambdaEI()->GetValue_f2()););

    // runaway/...
    DEF_RE_F1("runaway/nu_s_f1", "Slowing down frequency (on p1 flux grid) [s^-1]", qd->Store(nr_re,   (n1_re+1)*n2_re, this->cqtyRunaway->GetNuS()->GetValue_f1()););
    DEF_RE_F2("runaway/nu_s_f2", "Slowing down frequency (on p2 flux grid) [s^-1]", qd->Store(nr_re,   n1_re*(n2_re+1), this->cqtyRunaway->GetNuS()->GetValue_f2()););
    DEF_RE_F1("runaway/nu_D_f1", "Pitch-angle scattering frequency (on p1 flux grid) [s^-1]", qd->Store(nr_re,   (n1_re+1)*n2_re, this->cqtyRunaway->GetNuD()->GetValue_f1()););
    DEF_RE_F2("runaway/nu_D_f2", "Pitch-angle scattering frequency (on p2 flux grid) [s^-1]", qd->Store(nr_re,   n1_re*(n2_re+1), this->cqtyRunaway->GetNuD()->GetValue_f2()););
    DEF_RE_F1("runaway/lnLambda_ee_f1", "Coulomb logarithm for e-e collisions (on p1 flux grid)", qd->Store(nr_re,   (n1_re+1)*n2_re, this->cqtyRunaway->GetLnLambdaEE()->GetValue_f1()););
    DEF_RE_F2("runaway/lnLambda_ee_f2", "Coulomb logarithm for e-e collisions (on p2 flux grid)", qd->Store(nr_re,   n1_re*(n2_re+1), this->cqtyRunaway->GetLnLambdaEE()->GetValue_f2()););
    DEF_RE_F1("runaway/lnLambda_ei_f1", "Coulomb logarithm for e-i collisions (on p1 flux grid)", qd->Store(nr_re,   (n1_re+1)*n2_re, this->cqtyRunaway->GetLnLambdaEI()->GetValue_f1()););
    DEF_RE_F2("runaway/lnLambda_ei_f2", "Coulomb logarithm for e-i collisions (on p2 flux grid)", qd->Store(nr_re,   n1_re*(n2_re+1), this->cqtyRunaway->GetLnLambdaEI()->GetValue_f2()););


    // Declare groups of parameters (for registering
    // multiple parameters in one go)

    // Automatically add elements to the "fluid",
    // "hottail" and "runaway" groups
    for (auto qty : all_quantities) {
        if (qty->GetName().substr(0, 5) == "fluid")
            this->groups["fluid"].push_back(qty->GetName());
        else if (qty->GetName().substr(0, 7) == "hottail")
            this->groups["hottail"].push_back(qty->GetName());
        else if (qty->GetName().substr(0, 7) == "runaway")
            this->groups["runaway"].push_back(qty->GetName());
    }
    
    this->groups["nu_s"] = {
        "hottail/nu_s_f1", "hottail/nu_s_f2",
        "runaway/nu_s_f1", "runaway/nu_s_f2"
    };
    this->groups["nu_D"] = {
        "hottail/nu_D_f1", "hottail/nu_D_f2",
        "runaway/nu_D_f1", "runaway/nu_D_f2"
    };
    this->groups["lnLambda"] = {
        "hottail/lnLambda_ee_f1", "hottail/lnLambda_ee_f2",
        "hottail/lnLambda_ei_f1", "hottail/lnLambda_ei_f2",
        "runaway/lnLambda_ee_f1", "runaway/lnLambda_ee_f2",
        "runaway/lnLambda_ei_f1", "runaway/lnLambda_ei_f2"
    };
}

