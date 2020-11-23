/**
 * Load settings for the equation system preconditioner.
 */

#include "DREAM/DiagonalPreconditioner.hpp"
#include "DREAM/IO.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"


using namespace DREAM;
using namespace std;

#define MODULENAME "solver/preconditioner"

/**
 * Define options for the preconditioner.
 *
 * s: Settings object to define options in.
 */
void SimulationGenerator::DefinePreconditionerSettings(Settings *s) {
    s->DefineSetting(MODULENAME "/enabled", "Enable physics-based preconditioning", (bool)true);
    s->DefineSetting(MODULENAME "/names", "Names of unknowns to override scales for.", (const string)"");
    s->DefineSetting(MODULENAME "/equation_scales", "List of equation scales to use.", 0, (real_t*)nullptr);
    s->DefineSetting(MODULENAME "/unknown_scales", "List of unknown scales to use.", 0, (real_t*)nullptr);
}

/**
 * Load preconditioner settings and construct
 * a DiagoanlPreconditioner object.
 *
 * s:           Settings object to load settings from.
 * uqh:         Unknown quantity handler.
 * nontrivials: List of IDs of the non-trivial unknowns of the equation system.
 */
DiagonalPreconditioner *SimulationGenerator::LoadPreconditionerSettings(
    Settings *s, FVM::UnknownQuantityHandler *uqh,
    const vector<len_t>& nontrivials
) {
    bool enabled = s->GetBool(MODULENAME "/enabled");

    printf("Enabled: %s\n", (enabled ? "yes":"no"));

    if (!enabled)
        return nullptr;

    len_t neqn, nuqn;

    vector<string> uqtyNames = s->GetStringList(MODULENAME "/names");
    const real_t *eqn_scales = s->GetRealArray(MODULENAME "/equation_scales", 1, &neqn);
    const real_t *uqn_scales = s->GetRealArray(MODULENAME "/unknown_scales", 1, &nuqn);

    DiagonalPreconditioner *dp = new DiagonalPreconditioner(
        uqh, nontrivials
    );

    if (neqn != uqtyNames.size())
        throw SettingsException("DiagonalPreconditioner: Length of 'names' and 'equation_scales' do not match.");
    if (nuqn != uqtyNames.size())
        throw SettingsException("DiagonalPreconditioner: Length of 'names' and 'unknown_scales' do not match.");

    // Set scales...
    for (len_t i = 0; i < uqtyNames.size(); i++) {
        // Ignore non-existent unknowns...
        if (!uqh->HasUnknown(uqtyNames[i])) {
            DREAM::IO::PrintWarning(
                "%s: No unknown quantity with name '%s' in equation system. "
                "Ignoring tolerance settings...",
                MODULENAME, uqtyNames[i].c_str()
            );
            continue;
        }

        len_t id = uqh->GetUnknownID(uqtyNames[i]);

        dp->SetEquationScale(id, eqn_scales[i]);
        dp->SetUnknownScale(id, uqn_scales[i]);
    }

    return dp;
}

