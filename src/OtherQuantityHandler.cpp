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
#include "DREAM/Equations/Scalar/WallCurrentTerms.hpp"
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
    PostProcessor *postProcessor, RunawayFluid *REFluid, FVM::UnknownQuantityHandler *unknowns,
    std::vector<UnknownQuantityEquation*> *unknown_equations, IonHandler *ions,
    FVM::Grid *fluidGrid, FVM::Grid *hottailGrid, FVM::Grid *runawayGrid,
    FVM::Grid *scalarGrid, struct eqn_terms *oqty_terms
) : cqtyHottail(cqtyHottail), cqtyRunaway(cqtyRunaway),
    postProcessor(postProcessor), REFluid(REFluid), unknowns(unknowns), unknown_equations(unknown_equations),
    ions(ions), fluidGrid(fluidGrid), hottailGrid(hottailGrid), runawayGrid(runawayGrid),
    scalarGrid(scalarGrid), tracked_terms(oqty_terms) {

    id_Eterm = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_n_re  = unknowns->GetUnknownID(OptionConstants::UQTY_N_RE);
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_Wcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_jtot  = unknowns->GetUnknownID(OptionConstants::UQTY_J_TOT);
    id_psip  = unknowns->GetUnknownID(OptionConstants::UQTY_POL_FLUX);
    id_Ip    = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    id_psi_edge = unknowns->GetUnknownID(OptionConstants::UQTY_PSI_EDGE);
    id_psi_wall = unknowns->GetUnknownID(OptionConstants::UQTY_PSI_WALL);

    if (hottailGrid != nullptr) 
        id_f_hot = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);
    if (runawayGrid != nullptr) 
        id_f_re = unknowns->GetUnknownID(OptionConstants::UQTY_F_RE);

    this->DefineQuantities();
}


/**
 * Destructor.
 */
OtherQuantityHandler::~OtherQuantityHandler() {
    for (auto it = this->all_quantities.begin(); it != this->all_quantities.end(); it++)
        delete *it;
    
    delete this->tracked_terms;

    delete [] kineticVectorHot;
    delete [] kineticVectorRE;
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
            this->RegisterQuantity(*it, true);
        return true;
    } else return false;
}

/**
 * Register a new "other" quantity to keep track of.
 *
 * name:       Name of quantity to register.
 * ignorefail: If true, silently skips unrecognized other quantities.
 */
void OtherQuantityHandler::RegisterQuantity(const std::string& name, bool ignorefail) {
    OtherQuantity *oq = GetByName(name);

    if (oq == nullptr) {
        // Is the given name a group of parameters?
        // Try to register it!
        if (!ignorefail && !RegisterGroup(name))
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
    bool scalarCreated = false;

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
        } else if (!scalarCreated && oq->GetName().substr(0, 7) == "scalar/") {
            sf->CreateStruct(group+"scalar");
            scalarCreated = true;
        }
        oq->SaveSFile(sf, path);
    }
}

/********************************
 * IMPLEMENTATION OF QUANTITIES *
 ********************************/
/**
 * Define all other quantities.
 *
 * NOTE: This method is parsed by Sphinx when generating the online
 * documentation for DREAM. To ensure that all other quantities are accurately
 * identified and documented, the following rules should be followed:
 *
 * - The name and description of the quantity must be given on the same line
 * - The line defining the quantity must have "DEF_" as its first non-whitespace
 *   characters.
 *
 * (The code auto-generating the list of other quantities iterates through the
 * code line-by-line. If the first non-whitespace characters on the line are
 * "DEF_", then it recognizes the definition of an other quantity. Next, the
 * first two strings on the line are parsed (i.e. contents within double
 * quotation marks, "") and identified as (1) the quantity name and (2)
 * description)
 */
void OtherQuantityHandler::DefineQuantities() {
    // XXX here we assume that all momentum grids are the same
    const len_t nr_ht = (this->hottailGrid==nullptr ? 0 : this->hottailGrid->GetNr());
    const len_t n1_ht = (this->hottailGrid==nullptr ? 0 : this->hottailGrid->GetMomentumGrid(0)->GetNp1());
    const len_t n2_ht = (this->hottailGrid==nullptr ? 0 : this->hottailGrid->GetMomentumGrid(0)->GetNp2());

    const len_t nr_re = (this->runawayGrid==nullptr ? 0 : this->runawayGrid->GetNr());
    const len_t n1_re = (this->runawayGrid==nullptr ? 0 : this->runawayGrid->GetMomentumGrid(0)->GetNp1());
    const len_t n2_re = (this->runawayGrid==nullptr ? 0 : this->runawayGrid->GetMomentumGrid(0)->GetNp2());

    if(hottailGrid != nullptr)
        kineticVectorHot = new real_t[hottailGrid->GetNCells()];
    if(runawayGrid != nullptr)
        kineticVectorRE  = new real_t[runawayGrid->GetNCells()];

    // HELPER MACROS (to make definitions more compact)
    // Define on scalar grid
    #define DEF_SC(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), scalarGrid, 1, FVM::FLUXGRIDTYPE_DISTRIBUTION, [this](QuantityData *qd) {FUNC}));

    // Define on fluid grid
    #define DEF_FL(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), fluidGrid, 1, FVM::FLUXGRIDTYPE_DISTRIBUTION, [this](QuantityData *qd) {FUNC}));
    #define DEF_FL_FR(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), fluidGrid, 1, FVM::FLUXGRIDTYPE_RADIAL, [this](QuantityData *qd) {FUNC}));
    #define DEF_FL_MUL(NAME, MUL, DESC, FUNC) \
        this->all_quantities.push_back(new OtherQuantity((NAME), (DESC), fluidGrid, (MUL), FVM::FLUXGRIDTYPE_DISTRIBUTION, [this](QuantityData *qd) {FUNC}));

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
    DEF_FL("fluid/conductivity", "Electric conductivity in SI, Sauter formula (based on Braams)", qd->Store(this->REFluid->GetElectricConductivity()););
    DEF_FL("fluid/Eceff", "Effective critical electric field [V/m]", qd->Store(this->REFluid->GetEffectiveCriticalField()););
    DEF_FL("fluid/Ecfree", "Connor-Hastie threshold field (calculated with n=n_free) [V/m]", qd->Store(this->REFluid->GetConnorHastieField_COMPLETESCREENING()););
    DEF_FL("fluid/Ectot", "Connor-Hastie threshold field (calculated with n=n_tot) [V/m]", qd->Store(this->REFluid->GetConnorHastieField_NOSCREENING()););
    DEF_FL("fluid/EDreic", "Dreicer electric field [V/m]", qd->Store(this->REFluid->GetDreicerElectricField()););
    DEF_FL("fluid/GammaAva", "Avalanche growth rate [s^-1]", qd->Store(this->REFluid->GetAvalancheGrowthRate()););
    DEF_FL("fluid/gammaDreicer", "Dreicer runaway rate [s^-1 m^-3]", qd->Store(this->REFluid->GetDreicerRunawayRate()););
    DEF_FL("fluid/gammaCompton", "Compton runaway rate [s^-1 m^-3]", qd->Store(this->REFluid->GetComptonRunawayRate()););
    if(tracked_terms->n_re_hottail_rate != nullptr){
        DEF_FL("fluid/gammaHottail", "Hottail runaway rate [s^-1 m^-3]", qd->Store(tracked_terms->n_re_hottail_rate->GetRunawayRate()););
    }
    DEF_FL("fluid/gammaTritium", "Tritium runaway rate [s^-1 m^-3]", 
        const real_t *gt = this->REFluid->GetTritiumRunawayRate();
        real_t *v = qd->StoreEmpty();
        for (len_t ir = 0; ir < this->fluidGrid->GetNr(); ir++)
            v[ir] = gt[ir] * this->ions->GetTritiumDensity(ir);
    );

    // Magnetic ripple resonant momentum
    if (tracked_terms->f_hot_ripple_Dxx != nullptr) {
        len_t nModes = tracked_terms->f_hot_ripple_Dxx->GetNumberOfModes();
        DEF_FL_MUL("fluid/f_hot_ripple_pmn", nModes, "Magnetic ripple resonant momentum for f_hot [mc]", qd->Store(this->tracked_terms->f_hot_ripple_Dxx->GetResonantMomentum()[0]););
    }
    if (tracked_terms->f_re_ripple_Dxx != nullptr) {
        len_t nModes = tracked_terms->f_re_ripple_Dxx->GetNumberOfModes();
        DEF_FL_MUL("fluid/f_re_ripple_pmn", nModes, "Magnetic ripple resonant momentum for f_re [mc]", qd->Store(this->tracked_terms->f_re_ripple_Dxx->GetResonantMomentum()[0]););
    }
    // TODO at some point in the future, the mode numbers should be stored in
    // the RadialGrid instead of the RipplePitchScattering term, and so this
    // if statement would become unnecessary...
    if (tracked_terms->f_hot_ripple_Dxx != nullptr) {
        DEF_FL("fluid/ripple_m", "Magnetic ripple poloidal mode number", qd->Store(this->tracked_terms->f_hot_ripple_Dxx->GetPoloidalModeNumbers()););
        DEF_FL("fluid/ripple_n", "Magnetic ripple toroidal mode number", qd->Store(this->tracked_terms->f_hot_ripple_Dxx->GetToroidalModeNumbers()););
    } else if (tracked_terms->f_re_ripple_Dxx != nullptr) {
        DEF_FL("fluid/ripple_m", "Magnetic ripple poloidal mode number", qd->Store(this->tracked_terms->f_hot_ripple_Dxx->GetPoloidalModeNumbers()););
        DEF_FL("fluid/ripple_n", "Magnetic ripple toroidal mode number", qd->Store(this->tracked_terms->f_hot_ripple_Dxx->GetToroidalModeNumbers()););
    }
    DEF_FL("fluid/lnLambdaC", "Coulomb logarithm (relativistic)", qd->Store(this->REFluid->GetLnLambda()->GetLnLambdaC()););
    DEF_FL("fluid/lnLambdaT", "Coulomb logarithm (thermal)", qd->Store(this->REFluid->GetLnLambda()->GetLnLambdaT()););
    DEF_FL("fluid/pCrit", "Critical momentum for avalanche, compton and tritium (in units of mc)", qd->Store(this->REFluid->GetEffectiveCriticalRunawayMomentum()););
    if(tracked_terms->n_re_hottail_rate != nullptr){
        DEF_FL("fluid/pCritHottail", "Critical momentum for hottail (in units of mc)", qd->Store(tracked_terms->n_re_hottail_rate->GetHottailCriticalMomentum()););
    }
    DEF_FL("fluid/runawayRate", "Total runaway rate, dn_RE / dt [s^-1 m^-3]", qd->Store(this->postProcessor->GetRunawayRate()););
    DEF_FL("fluid/qR0", "Safety factor multiplied by major radius R0 [m]",
        real_t *vec = qd->StoreEmpty();
        const real_t *jtot = this->unknowns->GetUnknownData(id_jtot);
        for(len_t ir=0; ir<this->fluidGrid->GetNr(); ir++){
            real_t mu0Ip = Constants::mu0 * TotalPlasmaCurrentFromJTot::EvaluateIpInsideR(ir,this->fluidGrid->GetRadialGrid(),jtot);
            vec[ir] = this->fluidGrid->GetRadialGrid()->SafetyFactorNormalized(ir,mu0Ip);
        }
    )
    DEF_FL("fluid/tauEERel", "Relativistic electron collision time (4*pi*lnL*n_cold*r^2*c)^-1 [s]", qd->Store(this->REFluid->GetElectronCollisionTimeRelativistic()););
    DEF_FL("fluid/tauEETh", "Thermal electron collision time (tauEERel * [2T/mc^2]^1.5) [s]", qd->Store(this->REFluid->GetElectronCollisionTimeThermal()););
    
    // Power terms in heat equation
    if (tracked_terms->T_cold_ohmic != nullptr)
        DEF_FL("fluid/Tcold_ohmic", "Ohmic heating power density [J s^-1 m^-3]",
            real_t *Eterm = this->unknowns->GetUnknownData(this->id_Eterm);
            real_t *vec = qd->StoreEmpty();
            for(len_t ir=0; ir<this->fluidGrid->GetNr(); ir++)
                vec[ir] = 0;
            this->tracked_terms->T_cold_ohmic->SetVectorElements(vec, Eterm);
        );
    if (tracked_terms->T_cold_fhot_coll != nullptr)
        DEF_FL("fluid/Tcold_fhot_coll", "Collisional heating power density by f_hot [J s^-1 m^-3]",
            real_t *fhot = this->unknowns->GetUnknownData(id_f_hot);
            real_t *vec = qd->StoreEmpty();
            for(len_t ir=0; ir<this->fluidGrid->GetNr(); ir++)
                vec[ir] = 0;
            this->tracked_terms->T_cold_fhot_coll->SetVectorElements(vec, fhot);
        );
    if (tracked_terms->T_cold_fre_coll != nullptr)
        DEF_FL("fluid/Tcold_fre_coll", "Collisional heating power density by f_re [J s^-1 m^-3]",
            real_t *fre = this->unknowns->GetUnknownData(id_f_re);
            real_t *vec = qd->StoreEmpty();
            for(len_t ir=0; ir<this->fluidGrid->GetNr(); ir++)
                vec[ir] = 0;
            this->tracked_terms->T_cold_fre_coll->SetVectorElements(vec, fre);
        );
    if (tracked_terms->T_cold_nre_coll != nullptr)
        DEF_FL("fluid/Tcold_nre_coll", "Collisional heating power density by n_re [J s^-1 m^-3]",
            real_t *nre = this->unknowns->GetUnknownData(id_n_re);
            real_t *vec = qd->StoreEmpty();
            for(len_t ir=0; ir<this->fluidGrid->GetNr(); ir++)
                vec[ir] = 0;
            this->tracked_terms->T_cold_nre_coll->SetVectorElements(vec, nre);
        );
    
    if (tracked_terms->T_cold_transport != nullptr)
        DEF_FL("fluid/Tcold_transport", "Transported power density [J s^-1 m^-3]",
            real_t *Tcold = this->unknowns->GetUnknownData(this->id_Tcold);
            real_t *vec = qd->StoreEmpty();
            for(len_t ir=0; ir<this->fluidGrid->GetNr(); ir++)
                vec[ir] = 0;
            this->tracked_terms->T_cold_transport->SetVectorElements(vec, Tcold);
        );
    if (tracked_terms->T_cold_radiation != nullptr)
        DEF_FL("fluid/Tcold_radiation", "Radiated power density [J s^-1 m^-3]",
            real_t *ncold = this->unknowns->GetUnknownData(this->id_ncold);
            real_t *vec = qd->StoreEmpty();
            for(len_t ir=0; ir<this->fluidGrid->GetNr(); ir++)
                vec[ir] = 0;
            this->tracked_terms->T_cold_radiation->SetVectorElements(vec, ncold);
        );
    if (tracked_terms->T_cold_ion_coll != nullptr)
        DEF_FL("fluid/Tcold_ion_coll", "Collisional heating power density by ions [J s^-1 m^-3]",
            real_t *vec = qd->StoreEmpty();
            for(len_t ir=0; ir<this->fluidGrid->GetNr(); ir++)
                vec[ir] = 0;
            this->tracked_terms->T_cold_ion_coll->SetVectorElements(vec, nullptr);
        );

    if (tracked_terms->T_cold_transport) {
        DEF_FL_FR("fluid/Wcold_Tcold_Ar", "Net radial heat advection [m/s]",
            const real_t *Ar = this->unknown_equations->at(this->id_Wcold)->GetOperator(this->id_Tcold)->GetAdvectionCoeffR(0);
            qd->Store(Ar);
        );
        DEF_FL_FR("fluid/Wcold_Tcold_Drr", "Net radial heat diffusion [m/s]",
            const real_t *Drr = this->unknown_equations->at(this->id_Wcold)->GetOperator(this->id_Tcold)->GetDiffusionCoeffRR(0);
            qd->Store(Drr);
        );
    }

    /* TODO: come up with a condition to activate this term; for now it is inpractically expensive to evaluate
    DEF_FL("fluid/Tcold_radiationFromNuS", "Radiated power density predicted by the Hesslow screened nuS model [J s^-1 m^-3]",
        SlowingDownFrequency *nuS = this->REFluid->GetNuS();
        CollisionQuantity::collqty_settings settings_free;
        CollisionQuantity::collqty_settings settings_screened;
        settings_free.collfreq_type = OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_COMPLETELY_SCREENED;
        settings_screened.collfreq_type = OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED;
        settings_free.collfreq_mode = settings_screened.collfreq_mode = OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL;
        settings_free.lnL_type = settings_screened.lnL_type = OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;
        settings_free.bremsstrahlung_mode = settings_screened.bremsstrahlung_mode = OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_STOPPING_POWER;
        settings_free.screened_diffusion = settings_screened.screened_diffusion = OptionConstants::COLLQTY_SCREENED_DIFFUSION_MODE_MAXWELLIAN;


        std::function<real_t(len_t,real_t)> weightFunc = ([nuS, &settings_free, &settings_screened](len_t ir, real_t p)
        {
            real_t v = Constants::c * p/sqrt(1+p*p);
            return Constants::me*Constants::c*v*p*(nuS->evaluateAtP(ir,p,&settings_screened) - nuS->evaluateAtP(ir,p,&settings_free));
        });
        real_t *ncold = this->unknowns->GetUnknownData(this->id_ncold);
        real_t *Tcold = this->unknowns->GetUnknownData(this->id_Tcold);
        real_t *vec = qd->StoreEmpty();
        for(len_t ir=0; ir<this->fluidGrid->GetNr(); ir++)
            vec[ir] = integrateWeightedMaxwellian(ir, ncold[ir], Tcold[ir], weightFunc);
    );
    */

    DEF_FL("fluid/W_hot", "Energy density in f_hot [J m^-3]",
        real_t *vec = qd->StoreEmpty();
        if(hottailGrid != nullptr){
            const real_t *f_hot = this->unknowns->GetUnknownData(id_f_hot);
            const len_t nr = this->hottailGrid->GetNr();
            const real_t pThreshold = postProcessor->GetPThreshold();
            bool hasThreshold = (pThreshold != 0);
            const FVM::MomentQuantity::pThresholdMode pMode = postProcessor->GetPThresholdMode();
            
            len_t offset = 0;
            for(len_t ir=0; ir<nr; ir++){
                FVM::MomentumGrid *mg = this->hottailGrid->GetMomentumGrid(ir);
                const len_t n1 = mg->GetNp1();
                const len_t n2 = mg->GetNp2();
                for(len_t i=0; i<n1; i++){
                    real_t envelope = 1;
                    if(hasThreshold) 
                        envelope = FVM::MomentQuantity::ThresholdEnvelope(i, pThreshold, pMode, mg, unknowns->GetUnknownData(id_Tcold)[ir]);
                    for(len_t j=0; j<n2; j++){
                        real_t kineticEnergy = Constants::me * Constants::c * Constants::c * (mg->GetGamma(i,j)-1);
                        this->kineticVectorHot[offset + n1*j + i] = envelope * kineticEnergy * f_hot[offset + n1*j + i];
                    }
                }
                vec[ir] = this->hottailGrid->IntegralMomentumAtRadius(ir, this->kineticVectorHot+offset);
                offset += n1*n2;
            }
        }
    );
    DEF_FL("fluid/W_re", "Energy density in f_re [J m^-3]",
        real_t *vec = qd->StoreEmpty();
        if(runawayGrid != nullptr){
            const real_t *f_re = this->unknowns->GetUnknownData(id_f_re);
            const len_t nr = this->runawayGrid->GetNr();
            len_t offset = 0;
            for(len_t ir=0; ir<nr; ir++){
                FVM::MomentumGrid *mg = this->runawayGrid->GetMomentumGrid(ir);
                const len_t n1 = mg->GetNp1();
                const len_t n2 = mg->GetNp2();
                for(len_t i=0; i<n1; i++)
                    for(len_t j=0; j<n2; j++){
                        real_t kineticEnergy = Constants::me * Constants::c * Constants::c * (mg->GetGamma(i,j)-1);
                        this->kineticVectorRE[offset + n1*j + i] = kineticEnergy * f_re[offset + n1*j + i];
                    }
                vec[ir] = this->runawayGrid->IntegralMomentumAtRadius(ir, this->kineticVectorRE+offset);
                offset += n1*n2;
            }
        }
    );
    DEF_FL("fluid/Zeff", "Effective charge", qd->Store(this->REFluid->GetIonHandler()->GetZeff()););

    // hottail/...
    DEF_HT_F1("hottail/Ar", "Net radial advection on hot electron grid [m/s]",
        const real_t *const* Ar = this->unknown_equations->at(this->id_f_hot)->GetOperator(this->id_f_hot)->GetAdvectionCoeffR();
        qd->Store(nr_ht+1, n1_ht*n2_ht, Ar);
    );
    DEF_HT_F1("hottail/Ap1", "Net first momentum advection on hot electron grid [m/s]",
        const real_t *const* Ap = this->unknown_equations->at(this->id_f_hot)->GetOperator(this->id_f_hot)->GetAdvectionCoeff1();
        qd->Store(nr_ht, (n1_ht+1)*n2_ht, Ap);
    );
    DEF_HT_F1("hottail/Ap2", "Net second momentum advection on hot electron grid [m/s]",
        const real_t *const* Axi = this->unknown_equations->at(this->id_f_hot)->GetOperator(this->id_f_hot)->GetAdvectionCoeff1();
        qd->Store(nr_ht, n1_ht*(n2_ht+1), Axi);
    );
    DEF_HT_FR("hottail/Drr", "Net radial diffusion on hot electron grid [m/s]",
        const real_t *const* Drr = this->unknown_equations->at(this->id_f_hot)->GetOperator(this->id_f_hot)->GetDiffusionCoeffRR();
        qd->Store(nr_ht+1, n1_ht*n2_ht, Drr);
    );
    DEF_HT_F1("hottail/Dpp", "Net momentum-momentum diffusion on hot electron grid [m/s]",
        const real_t *const* Dpp = this->unknown_equations->at(this->id_f_hot)->GetOperator(this->id_f_hot)->GetDiffusionCoeff11();
        qd->Store(nr_ht, (n1_ht+1)*n2_ht, Dpp);
    );
    DEF_HT_F1("hottail/Dpx", "Net momentum-pitch diffusion on hot electron grid [m/s]",
        const real_t *const* Dpx = this->unknown_equations->at(this->id_f_hot)->GetOperator(this->id_f_hot)->GetDiffusionCoeff12();
        qd->Store(nr_ht, (n1_ht+1)*n2_ht, Dpx);
    );
    DEF_HT_F2("hottail/Dxp", "Net pitch-momentum diffusion on hot electron grid [m/s]",
        const real_t *const* Dxp = this->unknown_equations->at(this->id_f_hot)->GetOperator(this->id_f_hot)->GetDiffusionCoeff21();
        qd->Store(nr_ht, n1_ht*(n2_ht+1), Dxp);
    );
    DEF_HT_F2("hottail/Dxx", "Net pitch-pitch diffusion on hot electron grid [m/s]",
        const real_t *const* Dxx = this->unknown_equations->at(this->id_f_hot)->GetOperator(this->id_f_hot)->GetDiffusionCoeff22();
        qd->Store(nr_ht, n1_ht*(n2_ht+1), Dxx);
    );
    DEF_HT_F1("hottail/nu_s_f1", "Slowing down frequency (on p1 flux grid) [s^-1]", qd->Store(nr_ht,   (n1_ht+1)*n2_ht, this->cqtyHottail->GetNuS()->GetValue_f1()););
    DEF_HT_F2("hottail/nu_s_f2", "Slowing down frequency (on p2 flux grid) [s^-1]", qd->Store(nr_ht,   n1_ht*(n2_ht+1), this->cqtyHottail->GetNuS()->GetValue_f2()););
    DEF_HT_F1("hottail/nu_D_f1", "Pitch-angle scattering frequency (on p1 flux grid) [s^-1]", qd->Store(nr_ht,   (n1_ht+1)*n2_ht, this->cqtyHottail->GetNuD()->GetValue_f1()););
    DEF_HT_F2("hottail/nu_D_f2", "Pitch-angle scattering frequency (on p2 flux grid) [s^-1]", qd->Store(nr_ht,   n1_ht*(n2_ht+1), this->cqtyHottail->GetNuD()->GetValue_f2()););
    DEF_HT_F1("hottail/nu_par_f1", "Energy scattering frequency (on p1 flux grid) [s^-1]", qd->Store(nr_ht,   (n1_ht+1)*n2_ht, this->cqtyHottail->GetNuPar()->GetValue_f1()););
    DEF_HT_F2("hottail/nu_par_f2", "Energy scattering frequency (on p2 flux grid) [s^-1]", qd->Store(nr_ht,   n1_ht*(n2_ht+1), this->cqtyHottail->GetNuPar()->GetValue_f2()););
    DEF_HT_F1("hottail/lnLambda_ee_f1", "Coulomb logarithm for e-e collisions (on p1 flux grid)", qd->Store(nr_ht,   (n1_ht+1)*n2_ht, this->cqtyHottail->GetLnLambdaEE()->GetValue_f1()););
    DEF_HT_F2("hottail/lnLambda_ee_f2", "Coulomb logarithm for e-e collisions (on p2 flux grid)", qd->Store(nr_ht,   n1_ht*(n2_ht+1), this->cqtyHottail->GetLnLambdaEE()->GetValue_f2()););
    DEF_HT_F1("hottail/lnLambda_ei_f1", "Coulomb logarithm for e-i collisions (on p1 flux grid)", qd->Store(nr_ht,   (n1_ht+1)*n2_ht, this->cqtyHottail->GetLnLambdaEI()->GetValue_f1()););
    DEF_HT_F2("hottail/lnLambda_ei_f2", "Coulomb logarithm for e-i collisions (on p2 flux grid)", qd->Store(nr_ht,   n1_ht*(n2_ht+1), this->cqtyHottail->GetLnLambdaEI()->GetValue_f2()););

    // runaway/...
    DEF_RE_FR("runaway/Ar", "Net radial advection on runaway electron grid [m/s]",
        const real_t *const* Ar = this->unknown_equations->at(this->id_f_re)->GetOperator(this->id_f_re)->GetAdvectionCoeffR();
        qd->Store(nr_re+1, n1_re*n2_re, Ar);
    );
    DEF_RE_F1("runaway/Ap1", "Net first momentum advection on runaway electron grid [m/s]",
        const real_t *const* Ap = this->unknown_equations->at(this->id_f_re)->GetOperator(this->id_f_re)->GetAdvectionCoeff1();
        qd->Store(nr_re, (n1_re+1)*n2_re, Ap);
    );
    DEF_RE_F2("runaway/Ap2", "Net second momentum advection on runaway electron grid [m/s]",
        const real_t *const* Axi = this->unknown_equations->at(this->id_f_re)->GetOperator(this->id_f_re)->GetAdvectionCoeff2();
        qd->Store(nr_re, n1_re*(n2_re+1), Axi);
    );
    DEF_RE_FR("runaway/Drr", "Net radial diffusion on runaway electron grid [m/s]",
        const real_t *const* Drr = this->unknown_equations->at(this->id_f_re)->GetOperator(this->id_f_re)->GetDiffusionCoeffRR();
        qd->Store(nr_re+1, n1_re*n2_re, Drr);
    );
    DEF_RE_F1("runaway/Dpp", "Net momentum-momentum diffusion on runaway electron grid [m/s]",
        const real_t *const* Dpp = this->unknown_equations->at(this->id_f_re)->GetOperator(this->id_f_re)->GetDiffusionCoeff11();
        qd->Store(nr_re, (n1_re+1)*n2_re, Dpp);
    );
    DEF_RE_F1("runaway/Dpx", "Net momentum-pitch diffusion on runaway electron grid [m/s]",
        const real_t *const* Dpx = this->unknown_equations->at(this->id_f_re)->GetOperator(this->id_f_re)->GetDiffusionCoeff12();
        qd->Store(nr_re, (n1_re+1)*n2_re, Dpx);
    );
    DEF_RE_F2("runaway/Dxp", "Net pitch-momentum diffusion on runaway electron grid [m/s]",
        const real_t *const* Dxp = this->unknown_equations->at(this->id_f_re)->GetOperator(this->id_f_re)->GetDiffusionCoeff21();
        qd->Store(nr_re, n1_re*(n2_re+1), Dxp);
    );
    DEF_RE_F2("runaway/Dxx", "Net pitch-pitch diffusion on runaway electron grid [m/s]",
        const real_t *const* Dxx = this->unknown_equations->at(this->id_f_re)->GetOperator(this->id_f_re)->GetDiffusionCoeff22();
        qd->Store(nr_re, n1_re*(n2_re+1), Dxx);
    );
    DEF_RE_F1("runaway/nu_s_f1", "Slowing down frequency (on p1 flux grid) [s^-1]", qd->Store(nr_re,   (n1_re+1)*n2_re, this->cqtyRunaway->GetNuS()->GetValue_f1()););
    DEF_RE_F2("runaway/nu_s_f2", "Slowing down frequency (on p2 flux grid) [s^-1]", qd->Store(nr_re,   n1_re*(n2_re+1), this->cqtyRunaway->GetNuS()->GetValue_f2()););
    DEF_RE_F1("runaway/nu_D_f1", "Pitch-angle scattering frequency (on p1 flux grid) [s^-1]", qd->Store(nr_re,   (n1_re+1)*n2_re, this->cqtyRunaway->GetNuD()->GetValue_f1()););
    DEF_RE_F2("runaway/nu_D_f2", "Pitch-angle scattering frequency (on p2 flux grid) [s^-1]", qd->Store(nr_re,   n1_re*(n2_re+1), this->cqtyRunaway->GetNuD()->GetValue_f2()););
    DEF_RE_F1("runaway/nu_D_f1", "Energy scattering frequency (on p1 flux grid) [s^-1]", qd->Store(nr_re,   (n1_re+1)*n2_re, this->cqtyRunaway->GetNuPar()->GetValue_f1()););
    DEF_RE_F2("runaway/nu_D_f2", "Energy scattering frequency (on p2 flux grid) [s^-1]", qd->Store(nr_re,   n1_re*(n2_re+1), this->cqtyRunaway->GetNuPar()->GetValue_f2()););
    DEF_RE_F1("runaway/lnLambda_ee_f1", "Coulomb logarithm for e-e collisions (on p1 flux grid)", qd->Store(nr_re,   (n1_re+1)*n2_re, this->cqtyRunaway->GetLnLambdaEE()->GetValue_f1()););
    DEF_RE_F2("runaway/lnLambda_ee_f2", "Coulomb logarithm for e-e collisions (on p2 flux grid)", qd->Store(nr_re,   n1_re*(n2_re+1), this->cqtyRunaway->GetLnLambdaEE()->GetValue_f2()););
    DEF_RE_F1("runaway/lnLambda_ei_f1", "Coulomb logarithm for e-i collisions (on p1 flux grid)", qd->Store(nr_re,   (n1_re+1)*n2_re, this->cqtyRunaway->GetLnLambdaEI()->GetValue_f1()););
    DEF_RE_F2("runaway/lnLambda_ei_f2", "Coulomb logarithm for e-i collisions (on p2 flux grid)", qd->Store(nr_re,   n1_re*(n2_re+1), this->cqtyRunaway->GetLnLambdaEI()->GetValue_f2()););

    // scalar/..
    DEF_SC("scalar/radialloss_n_re", "Rate of runaway number loss through plasma edge, normalized to R0 [s^-1 m^-1]",
        const real_t *nre = this->unknowns->GetUnknownData(this->id_n_re);
        real_t v = 0;
        len_t nr = this->fluidGrid->GetNr();
        if (this->tracked_terms->n_re_advective_bc != nullptr)
            this->tracked_terms->n_re_advective_bc->AddToVectorElements((&v)-(nr-1), nre);
        if (this->tracked_terms->n_re_diffusive_bc != nullptr)
            this->tracked_terms->n_re_diffusive_bc->AddToVectorElements((&v)-(nr-1), nre);
        // multiply the flux through the boundary by the surface area (normalized to the major radius R0)
        v *= this->fluidGrid->GetVpVol(nr-1) * this->fluidGrid->GetRadialGrid()->GetDr(nr-1); 
        qd->Store(&v);
    );

    DEF_SC("scalar/energyloss_T_cold", "Rate of energy loss through plasma edge from T_cold transport, normalized to R0 [J s^-1 m^-1]",
        const real_t *Tcold = this->unknowns->GetUnknownData(this->id_Tcold);
        real_t v=0;
        len_t nr = this->fluidGrid->GetNr();
        if (this->tracked_terms->T_cold_advective_bc != nullptr)
            this->tracked_terms->T_cold_advective_bc->AddToVectorElements((&v)-(nr-1), Tcold);
        if (this->tracked_terms->T_cold_diffusive_bc != nullptr)
            this->tracked_terms->T_cold_diffusive_bc->AddToVectorElements((&v)-(nr-1), Tcold);
        
        // multiply the flux through the boundary by the surface area (normalized to the major radius R0)
        v *= this->fluidGrid->GetVpVol(nr-1) * this->fluidGrid->GetRadialGrid()->GetDr(nr-1);
        qd->Store(&v);
    );

    if(this->tracked_terms->f_re_advective_bc != nullptr || this->tracked_terms->f_re_diffusive_bc != nullptr ) {
        DEF_SC("scalar/radialloss_f_re", "Rate of particle number loss through plasma edge from f_re transport, normalized to R0 [s^-1 m^-1]",
            real_t v = integratedKineticBoundaryTerm( 
                this->id_f_re, [](len_t,len_t, FVM::MomentumGrid*){ return 1; }, 
                this->runawayGrid, this->tracked_terms->f_re_advective_bc, this->tracked_terms->f_re_diffusive_bc, kineticVectorRE 
            );
            qd->Store(&v);
        );
        DEF_SC("scalar/energyloss_f_re", "Rate of energy loss through plasma edge from f_re transport, normalized to R0 [J s^-1 m^-1]",
            real_t v = integratedKineticBoundaryTerm( 
                this->id_f_re, [](len_t i,len_t j, FVM::MomentumGrid *mg){ return Constants::me * Constants::c * Constants::c * (mg->GetGamma(i,j)-1); }, 
                this->runawayGrid, this->tracked_terms->f_re_advective_bc, this->tracked_terms->f_re_diffusive_bc, kineticVectorRE 
            );
            qd->Store(&v);
        );
    }
    if(this->tracked_terms->f_hot_advective_bc != nullptr || this->tracked_terms->f_hot_diffusive_bc != nullptr ) {
        DEF_SC("scalar/radialloss_f_hot", "Rate of particle number loss through plasma edge from f_hot transport, normalized to R0 [s^-1 m^-1]",
            real_t v = integratedKineticBoundaryTerm( 
                this->id_f_hot, [](len_t,len_t, FVM::MomentumGrid*){ return 1; }, 
                this->hottailGrid, this->tracked_terms->f_hot_advective_bc, this->tracked_terms->f_hot_diffusive_bc, kineticVectorHot
            );
            qd->Store(&v);
        );
        
        DEF_SC("scalar/energyloss_f_hot", "Rate of energy loss through plasma edge from f_hot transport, normalized to R0 [J s^-1 m^-1]",
            real_t v = integratedKineticBoundaryTerm( 
                this->id_f_hot, [](len_t i,len_t j, FVM::MomentumGrid *mg){ return Constants::me * Constants::c * Constants::c * (mg->GetGamma(i,j)-1); }, 
                this->hottailGrid, this->tracked_terms->f_hot_advective_bc, this->tracked_terms->f_hot_diffusive_bc, kineticVectorHot
            );
            qd->Store(&v);
        );
    }

    // Magnetic energy and internal inductance
    DEF_SC("scalar/E_mag", "Total energy contained in the poloidal magnetic field within the vessel, normalized to R0 [J/m]",
        real_t v = evaluateMagneticEnergy();
        qd->Store(&v);
    );
    DEF_SC("scalar/L_i", "Internal inductance for poloidal magnetic energy normalized to R0 [J/A^2 m]",
        const real_t Ip = this->unknowns->GetUnknownData(id_Ip)[0];
        real_t v = 2*evaluateMagneticEnergy() / (Ip*Ip);
        qd->Store(&v);
    );
    DEF_SC("scalar/l_i", "Normalized internal inductance for poloidal magnetic energy (2Li/mu0R0)",
        const real_t Ip = this->unknowns->GetUnknownData(id_Ip)[0];
        real_t Li = 2*evaluateMagneticEnergy() / (Ip*Ip);
        real_t v = Li * 2/Constants::mu0;
        qd->Store(&v);
    );
    DEF_SC("scalar/L_i_flux", "Internal inductance for poloidal flux psi_p, normalized to R0 [J/A^2 m]",
        const real_t Ip = this->unknowns->GetUnknownData(id_Ip)[0];
        const real_t psip_0 = this->unknowns->GetUnknownData(id_psip)[0];
        const real_t psip_a = this->unknowns->GetUnknownData(id_psi_edge)[0];
        real_t v = (psip_a - psip_0) / Ip;
        qd->Store(&v);
    );
    
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
        else if (qty->GetName().substr(0, 6) == "scalar")
            this->groups["scalar"].push_back(qty->GetName());
    }

    this->groups["ripple"] = {
        "fluid/ripple_m", "fluid/ripple_n", "fluid/f_hot_ripple_pmn", "fluid/f_re_ripple_pmn"
    };
    this->groups["transport"] = {
        "scalar/radialloss_n_re", "scalar/energyloss_T_cold", 
        "scalar/radialloss_f_re", "scalar/energyloss_f_re",
        "scalar/radialloss_f_hot", "scalar/energyloss_f_hot"
    };
    this->groups["nu_s"] = {
        "hottail/nu_s_f1", "hottail/nu_s_f2",
        "runaway/nu_s_f1", "runaway/nu_s_f2"
    };
    this->groups["nu_D"] = {
        "hottail/nu_D_f1", "hottail/nu_D_f2",
        "runaway/nu_D_f1", "runaway/nu_D_f2"
    };
    this->groups["nu_par"] = {
        "hottail/nu_par_f1", "hottail/nu_par_f2",
        "runaway/nu_par_f1", "runaway/nu_par_f2"
    };
    this->groups["lnLambda"] = {
        "hottail/lnLambda_ee_f1", "hottail/lnLambda_ee_f2",
        "hottail/lnLambda_ei_f1", "hottail/lnLambda_ei_f2",
        "runaway/lnLambda_ee_f1", "runaway/lnLambda_ee_f2",
        "runaway/lnLambda_ei_f1", "runaway/lnLambda_ei_f2"
    };
    this->groups["energy"] = {
        "fluid/Tcold_ohmic", "fluid/Tcold_fhot_coll", "fluid/Tcold_fre_coll",
        "fluid/Tcold_transport", "fluid/Tcold_radiation", "fluid/Tcold_ion_coll",
        "fluid/W_hot", "fluid/W_re", "scalar/E_mag", "scalar/L_i", "scalar/l_i",
        "scalar/energyloss_T_cold", "scalar/energyloss_f_re", "scalar/energyloss_f_hot"
    };
}


/**
 * Returns the scalar corresponding a kinetic boundary condition term which has 
 * been integrated over momentum with a provided weight function.
 * 
 * Parameters
 *  id_f:           unknown id of kinetic quantity
 *  momentFunction: weight function that we integrate the equation term over
 *  grid:           grid on which the kinetic quantity lives
 *  advective_bc:   the first boundary condition equation term
 *  diffusive_bc:   the second boundary condition equation term
 *  kineticVector:  an array sufficiently large to contain the kinetic grid 
 */
real_t OtherQuantityHandler::integratedKineticBoundaryTerm(
        len_t id_f, std::function<real_t(len_t,len_t,FVM::MomentumGrid*)> momentFunction, FVM::Grid *grid, 
        FVM::BC::BoundaryCondition *advective_bc, FVM::BC::BoundaryCondition *diffusive_bc, 
        real_t *kineticVector
) {
    const real_t *f = this->unknowns->GetUnknownData(id_f); 
    len_t nr = grid->GetNr();
    FVM::MomentumGrid *mg = grid->GetMomentumGrid(nr-1);
    len_t n1 = mg->GetNp1();
    len_t n2 = mg->GetNp2();
    len_t offset = grid->GetNCells() - n1*n2; // get offset for ir=nr-1
    for(len_t i=0; i<n1*n2; i++) // reset vector
        kineticVector[offset+i] = 0;
    if (advective_bc != nullptr)
        advective_bc->AddToVectorElements(kineticVector,f);
    if (diffusive_bc != nullptr)
        diffusive_bc->AddToVectorElements(kineticVector,f);

    // take energy moment of the boundary condition
    for(len_t i=0; i<n1; i++)
        for(len_t j=0; j<n2; j++)
            kineticVector[offset+n1*j+i] *= momentFunction(i,j,mg);

    real_t v = grid->IntegralMomentumAtRadius(nr-1,kineticVector+offset);
    v *= grid->GetVpVol(nr-1) * grid->GetRadialGrid()->GetDr(nr-1);
    return v;
}

/** 
 * Returns the total poloidal magnetic energy internal 
 * to the tokamak chamber normalized to R0
 */
real_t OtherQuantityHandler::evaluateMagneticEnergy(){
    FVM::RadialGrid *rGrid = this->fluidGrid->GetRadialGrid();
    const real_t *G_R0 = rGrid->GetBTorG();
    const real_t *VpVol = rGrid->GetVpVol();
    const real_t *dr = rGrid->GetDr();
    const real_t *FSA_1OverR2 = rGrid->GetFSA_1OverR2();
    const real_t *Bmin = rGrid->GetBmin();
    const real_t *jtot = this->unknowns->GetUnknownData(id_jtot);
    const real_t *psi_p = this->unknowns->GetUnknownData(id_psip);
    const real_t psi_p_wall = this->unknowns->GetUnknownData(id_psi_wall)[0];
    const real_t Ip = this->unknowns->GetUnknownData(id_Ip)[0];
    real_t E_mag = .5 * psi_p_wall*Ip;
    real_t fourPiInv = 1/(4*M_PI);
    for(len_t ir=0; ir<rGrid->GetNr(); ir++)
        E_mag -= fourPiInv*dr[ir] * VpVol[ir] * G_R0[ir] * FSA_1OverR2[ir] * jtot[ir] * psi_p[ir] / Bmin[ir];
    
    return E_mag;
}

/** 
 * GSL function definitions defining the integrand of the Maxwellian moment.
 * Used in 'OtherQuantityHandler::integrateWeightedMaxwellian'
 */
struct MaxwellianIntegrandParams {len_t ir; real_t n; real_t T; std::function<real_t(len_t,real_t)> weightFunc;};
real_t MaxwellianIntegrandFunc(real_t p, void *par){
    MaxwellianIntegrandParams *params = (MaxwellianIntegrandParams*) par;
    return 4*M_PI*p*p*params->weightFunc(params->ir, p)*Constants::RelativisticMaxwellian(p,params->n, params->T);
}

/**
 * Evaluates the 'WeightFunc' (angle-averaged) moment over 
 * a relativistic Maxwellian at density n and temperature T.
 * 'WeightFunc(ir,p)' is a function of radial grid index and momentum.
 * Integrates adaptively from 0 to infinity
 */
real_t OtherQuantityHandler::integrateWeightedMaxwellian(len_t ir, real_t n, real_t T, std::function<real_t(len_t,real_t)> weightFunc){
    gsl_integration_workspace *gsl_ad_w = gsl_integration_workspace_alloc(1000);
    
    MaxwellianIntegrandParams params = {ir,n,T,weightFunc};
    gsl_function GSL_Func;
    GSL_Func.function = &(MaxwellianIntegrandFunc);
    GSL_Func.params = &params;
    real_t result, error, reltol=1e-4, abstol=0;
    gsl_integration_qagiu(&GSL_Func, 0, abstol, reltol, gsl_ad_w->limit, gsl_ad_w, &result, &error);
    gsl_integration_workspace_free(gsl_ad_w);

    return result;
}
