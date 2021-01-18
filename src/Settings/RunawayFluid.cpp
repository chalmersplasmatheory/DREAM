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
RunawayFluid *SimulationGenerator::ConstructRunawayFluid(FVM::Grid *g,
    FVM::UnknownQuantityHandler *unknowns, IonHandler *ih, 
    AnalyticDistributionRE *distRE, AnalyticDistributionHottail *distHT,
    OptionConstants::momentumgrid_type gridtype, Settings *s
) {
    struct CollisionQuantity::collqty_settings *cq =
        new CollisionQuantity::collqty_settings;

    cq->collfreq_mode       = OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
    cq->collfreq_type       = (enum OptionConstants::collqty_collfreq_type)s->GetInteger(MODNAME "/collfreq_type");
    cq->lnL_type            = (enum OptionConstants::collqty_lnLambda_type)s->GetInteger(MODNAME "/lnlambda");
    cq->pstar_mode          = (enum OptionConstants::collqty_pstar_mode)s->GetInteger(MODNAME "/pstar_mode");
    cq->screened_diffusion  = (enum OptionConstants::collqty_screened_diffusion_mode)s->GetInteger(MODNAME "/screened_diffusion_mode");

    OptionConstants::conductivity_mode cond_mode      = (enum OptionConstants::conductivity_mode)    s->GetInteger("eqsys/j_ohm/conductivityMode");
    OptionConstants::eqterm_dreicer_mode dreicer_mode = (enum OptionConstants::eqterm_dreicer_mode)  s->GetInteger("eqsys/n_re/dreicer");
    OptionConstants::collqty_Eceff_mode Eceff_mode    = (enum OptionConstants::collqty_Eceff_mode)   s->GetInteger("eqsys/n_re/Eceff");
    OptionConstants::eqterm_avalanche_mode ava_mode   = (enum OptionConstants::eqterm_avalanche_mode)s->GetInteger("eqsys/n_re/avalanche");
    OptionConstants::eqterm_compton_mode compton_mode = (enum OptionConstants::eqterm_compton_mode)  s->GetInteger("eqsys/n_re/compton/mode");
    real_t compton_photon_flux = s->GetReal("eqsys/n_re/compton/flux");

    // Note: these collision quantities will only be used for their evaluateAt(..., inSettings) 
    //       methods inside REFluid, and be called with other settings than 'cq'. 
    CoulombLogarithm *lnLEE = new CoulombLogarithm(g,unknowns,ih,gridtype,cq,CollisionQuantity::LNLAMBDATYPE_EE);
    CoulombLogarithm *lnLEI = new CoulombLogarithm(g,unknowns,ih,gridtype,cq,CollisionQuantity::LNLAMBDATYPE_EI);
    SlowingDownFrequency *nuS  = new SlowingDownFrequency(g,unknowns,ih,lnLEE,lnLEI,gridtype,cq);
    PitchScatterFrequency *nuD = new PitchScatterFrequency(g,unknowns,ih,lnLEI,lnLEE,gridtype,cq);

    real_t thresholdToNeglectTrapped = 100*sqrt(std::numeric_limits<real_t>::epsilon());
    OptionConstants::eqterm_hottail_dist_mode ht_dist_mode = (enum OptionConstants::eqterm_hottail_dist_mode)s->GetInteger("eqsys/f_hot/hottailDist");
    distRE = new AnalyticDistributionRE(g->GetRadialGrid(), nuD, Eceff_mode, thresholdToNeglectTrapped);
    distHT = new AnalyticDistributionHottail(g->GetRadialGrid(), unknowns, ht_dist_mode);
    RunawayFluid *REF = new RunawayFluid(
        g, unknowns, nuS,nuD,lnLEE,lnLEI, cq, ih, distRE, cond_mode, dreicer_mode,
        Eceff_mode, ava_mode, compton_mode, compton_photon_flux
    );
    
    return REF;
}

