/**
 * Set up a collision quantity handler.
 */

#include "DREAM/Settings/SimulationGenerator.hpp"


using namespace DREAM;
using namespace std;


#define MODNAME "collisions"

/**
 * Construct a CollisionQuantityHandler object for the specified
 * grid object. Read settings from the named section of the settings.
 *
 * name:     Name of settings section to load options from.
 * grid:     Grid object for which to construct the collision handler.
 * unknowns: List of unknowns in the associated equation system.
 * s:        Settings describing how to construct the collision handler.
 */
void SimulationGenerator::ConstructRunawayFluid(FVM::Grid *g,
    FVM::UnknownQuantityHandler *unknowns, IonHandler *ih,
    OptionConstants::momentumgrid_type gridtype, 
    EquationSystem *eqsys, Settings *s
) {
    struct CollisionQuantity::collqty_settings *cqsetForPc = new CollisionQuantity::collqty_settings;
    struct CollisionQuantity::collqty_settings *cqsetForEc = new CollisionQuantity::collqty_settings;

    cqsetForPc->collfreq_mode      = cqsetForEc->collfreq_mode      = OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
    cqsetForPc->collfreq_type      = cqsetForEc->collfreq_type      = (enum OptionConstants::collqty_collfreq_type)          s->GetInteger(MODNAME "/collfreq_type");
    cqsetForPc->screened_diffusion = cqsetForEc->screened_diffusion = (enum OptionConstants::collqty_screened_diffusion_mode)s->GetInteger(MODNAME "/screened_diffusion_mode");
    cqsetForPc->pstar_mode         = cqsetForEc->pstar_mode         = (enum OptionConstants::collqty_pstar_mode)             s->GetInteger(MODNAME "/pstar_mode");

    cqsetForPc->lnL_type = (enum OptionConstants::collqty_lnLambda_type)s->GetInteger(MODNAME "/lnlambda");
    cqsetForEc->lnL_type = OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT; 

    cqsetForPc->bremsstrahlung_mode = OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT;
    cqsetForEc->bremsstrahlung_mode = OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_STOPPING_POWER;

    OptionConstants::conductivity_mode cond_mode      = (enum OptionConstants::conductivity_mode)    s->GetInteger("eqsys/j_ohm/conductivityMode");
    OptionConstants::eqterm_dreicer_mode dreicer_mode = (enum OptionConstants::eqterm_dreicer_mode)  s->GetInteger("eqsys/n_re/dreicer");
    OptionConstants::collqty_Eceff_mode Eceff_mode    = (enum OptionConstants::collqty_Eceff_mode)   s->GetInteger("eqsys/n_re/Eceff");
    OptionConstants::eqterm_avalanche_mode ava_mode   = (enum OptionConstants::eqterm_avalanche_mode)s->GetInteger("eqsys/n_re/avalanche");
    OptionConstants::eqterm_compton_mode compton_mode = (enum OptionConstants::eqterm_compton_mode)  s->GetInteger("eqsys/n_re/compton/mode");
    FVM::Interpolator1D *compton_photon_flux = LoadDataT("eqsys/n_re/compton", s, "flux");

    // Note: these collision quantities will only be used for their evaluateAt(..., inSettings) 
    //       methods inside REFluid, and be called with other settings than 'cq'. 
    CoulombLogarithm *lnLEE = new CoulombLogarithm(g,unknowns,ih,gridtype,cqsetForPc,CollisionQuantity::LNLAMBDATYPE_EE);
    CoulombLogarithm *lnLEI = new CoulombLogarithm(g,unknowns,ih,gridtype,cqsetForPc,CollisionQuantity::LNLAMBDATYPE_EI);
    SlowingDownFrequency *nuS  = new SlowingDownFrequency(g,unknowns,ih,lnLEE,lnLEI,gridtype,cqsetForPc);
    PitchScatterFrequency *nuD = new PitchScatterFrequency(g,unknowns,ih,lnLEI,lnLEE,gridtype,cqsetForPc);

    real_t thresholdToNeglectTrapped = 100*sqrt(std::numeric_limits<real_t>::epsilon());
    OptionConstants::uqty_f_hot_dist_mode ht_dist_mode = (enum OptionConstants::uqty_f_hot_dist_mode)s->GetInteger("eqsys/f_hot/dist_mode");
    AnalyticDistributionRE::dist_mode re_dist_mode = (Eceff_mode==OptionConstants::COLLQTY_ECEFF_MODE_SIMPLE) ? 
            AnalyticDistributionRE::RE_PITCH_DIST_SIMPLE : AnalyticDistributionRE::RE_PITCH_DIST_FULL;
    AnalyticDistributionRE *distRE = new AnalyticDistributionRE(g->GetRadialGrid(), unknowns, nuD, cqsetForEc, re_dist_mode, thresholdToNeglectTrapped);
    AnalyticDistributionHottail *distHT = nullptr;
    OptionConstants::uqty_distribution_mode fhot_mode = (enum OptionConstants::uqty_distribution_mode)s->GetInteger("eqsys/f_hot/mode");
    if(fhot_mode == OptionConstants::UQTY_DISTRIBUTION_MODE_ANALYTICAL){
        FVM::RadialGrid *rGrid = g->GetRadialGrid();
        real_t *n0 = LoadDataR("eqsys/f_hot", rGrid, s, "n0");
        real_t *T0 = LoadDataR("eqsys/f_hot", rGrid, s, "T0");
        distHT = new AnalyticDistributionHottail(rGrid, unknowns, n0, T0, ht_dist_mode);
    }
    
    RunawayFluid *REF = new RunawayFluid(
        g, unknowns, nuS, nuD, lnLEE, lnLEI, ih, distRE, cqsetForPc, cqsetForEc,
        cond_mode,dreicer_mode,Eceff_mode,ava_mode,compton_mode,compton_photon_flux
    );
    distRE->SetREFluid(REF);
    eqsys->SetAnalyticDists(distRE, distHT);
    eqsys->SetREFluid(REF);
}

