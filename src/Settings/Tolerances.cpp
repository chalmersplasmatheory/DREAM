/**
 * Load tolerance settings.
 */

#include "DREAM/IO.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"


using namespace DREAM;
using namespace std;


/**
 * Define options for tolerance settings.
 *
 * modname: Name of module to place the settings in.
 * s:       Settings object to define options in.
 * name:    Name of tolerance settings object to use.
 */
void SimulationGenerator::DefineToleranceSettings(
    const string& modname, Settings *s,
    const string& name
) {
    const string n = modname + "/" + name;

    s->DefineSetting(n + "/reltol",  "General relative tolerance to apply.", (real_t)1e-6);
    s->DefineSetting(n + "/names",   "Names of unknowns to override tolerances for", (const string)"");
    s->DefineSetting(n + "/abstols", "List of absolute tolerances to apply to the unknowns in 'names'.", 0, (real_t*)nullptr);
    s->DefineSetting(n + "/reltols", "List of relative tolerances to apply to the unknowns in 'names'.", 0, (real_t*)nullptr);
}

/**
 * Load tolerance settings for the specified module.
 * 
 * modname:     Name of module to place the settings in.
 * s:           Settings object to define options in.
 * uqh:         Unknown quantity handler.
 * nontrivials: List of IDs of the non-trivial unknowns of the equation system.
 * name:        Name of tolerance settings object to use.
 */
ConvergenceChecker *SimulationGenerator::LoadToleranceSettings(
    const string& modname, Settings *s,
	vector<UnknownQuantityEquation*> *unknown_equations,
    FVM::UnknownQuantityHandler *uqh, const vector<len_t>& nontrivials,
    const string& name
) {
    const string n = modname + "/" + name;

    len_t nabstols, nreltols;

    real_t reltol = s->GetReal(n + "/reltol");
    vector<string> uqtyNames = s->GetStringList(n + "/names");
    const real_t *abstols = s->GetRealArray(n + "/abstols", 1, &nabstols);
    const real_t *reltols = s->GetRealArray(n + "/reltols", 1, &nreltols);

    ConvergenceChecker *cc = new ConvergenceChecker(
        uqh, unknown_equations, nontrivials, nullptr, reltol
    );

    // Set tolerances...
    for (len_t i = 0; i < uqtyNames.size(); i++) {
        // Ignore non-existent unknowns...
        // (These settings may be remnants from previous settings
        //  objects, and so instead of stopping the simulation we
        //  simply ignore it and emit a warning in case it's a typo)
        if (!uqh->HasUnknown(uqtyNames[i])) {
            DREAM::IO::PrintWarning(
                "%s: No unknown quantity with name '%s' in equation system. "
                "Ignoring tolerance settings...",
                n.c_str(), uqtyNames[i].c_str()
            );
            continue;
        }

        len_t id = uqh->GetUnknownID(uqtyNames[i]);

        cc->SetAbsoluteTolerance(id, abstols[i]);
        cc->SetRelativeTolerance(id, reltols[i]);
    }

    return cc;
}

