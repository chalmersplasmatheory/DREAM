/**
 * Set up a collision quantity handler.
 */

#include <string>
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"


using namespace DREAM;
using namespace std;


#define MODNAME "collisions"

/**
 * Define options applying to collision models.
 *
 * name: Name of grid settings group to define options in.
 * s:    Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_CollisionQuantityHandler(
    Settings *s
) {
    /*s->DefineSetting(mod + "/" MODNAME "/lnlambda", "Model to use when evaluating Coulomb logarithm", (int_t)OptionConstants::COLLQTY_LNLAMBDA_CONSTANT);
    s->DefineSetting(mod + "/" MODNAME "/collfreq_mode", "Mode in which to evaluate collision frequencies", (int_t)OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL);
    s->DefineSetting(mod + "/" MODNAME "/collfreq_type", "Model to use when evaluating collision frequencies", (int_t)OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED);
    s->DefineSetting(mod + "/" MODNAME "/bremsstrahlung", "Model to use for bremsstrahlung", (int_t)OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT);*/
    s->DefineSetting(MODNAME "/lnlambda", "Model to use when evaluating Coulomb logarithm", (int_t)OptionConstants::COLLQTY_LNLAMBDA_CONSTANT);
    s->DefineSetting(MODNAME "/collfreq_mode", "Mode in which to evaluate collision frequencies", (int_t)OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL);
    s->DefineSetting(MODNAME "/collfreq_type", "Model to use when evaluating collision frequencies", (int_t)OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED);
    s->DefineSetting(MODNAME "/bremsstrahlung_mode", "Model to use for bremsstrahlung", (int_t)OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT);
    s->DefineSetting(MODNAME "/pstar_mode", "Model to use for p_\\star", (int_t)OptionConstants::COLLQTY_PSTAR_MODE_COLLISIONLESS);
    s->DefineSetting(MODNAME "/screened_diffusion_mode", "Model to use for the energy diffusion frequency caused by bound electrons", (int_t)OptionConstants::COLLQTY_SCREENED_DIFFUSION_MODE_MAXWELLIAN);
}

/**
 * Construct a CollisionQuantityHandler object for the specified
 * grid object. Read settings from the named section of the settings.
 *
 * name:     Name of settings section to load options from.
 * grid:     Grid object for which to construct the collision handler.
 * unknowns: List of unknowns in the associated equation system.
 * s:        Settings describing how to construct the collision handler.
 */
CollisionQuantityHandler *SimulationGenerator::ConstructCollisionQuantityHandler(
    enum OptionConstants::momentumgrid_type gridtype, FVM::Grid *grid,
    FVM::UnknownQuantityHandler *unknowns, IonHandler *ionHandler,  Settings *s
) {
    struct CollisionQuantity::collqty_settings *cq =
        new CollisionQuantity::collqty_settings;

    cq->collfreq_type       = (enum OptionConstants::collqty_collfreq_type)           s->GetInteger(MODNAME "/collfreq_type");
    cq->collfreq_mode       = (enum OptionConstants::collqty_collfreq_mode)           s->GetInteger(MODNAME "/collfreq_mode");
    cq->lnL_type            = (enum OptionConstants::collqty_lnLambda_type)           s->GetInteger(MODNAME "/lnlambda");
    cq->bremsstrahlung_mode = (enum OptionConstants::eqterm_bremsstrahlung_mode)      s->GetInteger(MODNAME "/bremsstrahlung_mode");
    cq->pstar_mode          = (enum OptionConstants::collqty_pstar_mode)              s->GetInteger(MODNAME "/pstar_mode");
    cq->screened_diffusion  = (enum OptionConstants::collqty_screened_diffusion_mode) s->GetInteger(MODNAME "/screened_diffusion_mode");

    CollisionQuantityHandler *cqh = new CollisionQuantityHandler(grid, unknowns, ionHandler,gridtype,cq);

    return cqh;
}

