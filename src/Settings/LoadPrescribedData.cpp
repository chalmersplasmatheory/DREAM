/**
 * Routines for loading prescribed data that should be put
 * on the computational grid (e.g. as a PrescribedParameter or
 * initial condition for some equation.
 *
 * All public routines in this file take a "module name" as input
 * and read data from the corresponding part of the Settings. If
 * the module name is "equation", then this module expects the
 * following settings to be present in the given 'Settings' object:
 *
 * RADIAL GRID
 *   equations/data/r    -- Radial grid on which the data is defined (nr elements)
 *   equations/data/x    -- Prescribed data (nr elements)
 *
 * RADIAL GRID + TIME GRID
 *   equations/data/r    -- Radial grid on which the data is defined (nr elements)
 *   equations/data/t    -- Time grid on which the data is defined (nt elements)
 *   equations/data/x    -- Prescribed data (size nt-by-nr)
 *
 * P/XI GRID
 *   equations/data/p    -- Momentum grid on which the data is defined (np elements)
 *   equations/data/x    -- Prescribed data (size nxi-by-np)
 *   equations/data/xi   -- Pitch grid on which the data is defined (nxi elements)
 */

#include <gsl/gsl_interp.h>
#include <string>
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"


using namespace DREAM;
using namespace std;


/**
 * Define options for a "radius+time" 'data' section in the
 * specified module.
 */
void SimulationGenerator::DefineDataRT(
    const string& modname, Settings *s
) {
    const len_t ndim[2] = {0,0};

    s->DefineSetting(modname + "/data/r", "Radial grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/data/rinterp", "Interpolation method to use for radial grid interpolation.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_GSL_LINEAR);
    s->DefineSetting(modname + "/data/t", "Time grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/data/tinterp", "Interpolation method to use for time grid interpolation.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_GSL_LINEAR);
    s->DefineSetting(modname + "/data/x", "Prescribed data.", 2, ndim, (real_t*)nullptr);
}

/**
 * Load data from the 'data' section of the specified module.
 * The data is expected to depend on radius and time. It is
 * interpolated in radius to the given radial grid.
 *
 * modname: Name of module to load data from.
 * rgrid:   Radial grid to interpolate data to.
 * s:       Settings object to load data from.
 */
FVM::Interpolator1D *SimulationGenerator::LoadDataRT(
    const string& modname, FVM::RadialGrid *rgrid, Settings *s
) {
    len_t xdims[2], nr_inp, nt;

    real_t *r = s->GetRealArray(modname + "/data/r", 1, &nr_inp);
    real_t *t = s->GetRealArray(modname + "/data/t", 1, &nt);
    real_t *x = s->GetRealArray(modname + "/data/x", 2, xdims);

    enum OptionConstants::prescribed_data_interp tinterp =
        (enum OptionConstants::prescribed_data_interp)s->GetInteger(modname + "/data/tinterp");
    enum OptionConstants::prescribed_data_interp_gsl rinterp =
        (enum OptionConstants::prescribed_data_interp_gsl)s->GetInteger(modname + "/data/rinterp");

    // Select Interpolator1D interpolation method
    enum FVM::Interpolator1D::interp_method interp1_meth;
    switch (tinterp) {
        case OptionConstants::PRESCRIBED_DATA_INTERP_NEAREST:
            interp1_meth = FVM::Interpolator1D::INTERP_NEAREST; break;
        case OptionConstants::PRESCRIBED_DATA_INTERP_LINEAR:
            interp1_meth = FVM::Interpolator1D::INTERP_LINEAR; break;

        default:
            throw SettingsException(
                "%s: Unrecognized interpolation method on time grid: %d.",
                modname.c_str(), tinterp
            );
    }

    // Select GSL interpolation method
    const gsl_interp_type *gsl_meth;
    switch (rinterp) {
        case OptionConstants::PRESCRIBED_DATA_INTERP_GSL_LINEAR:
            gsl_meth = gsl_interp_linear; break;
        case OptionConstants::PRESCRIBED_DATA_INTERP_GSL_POLYNOMIAL:
            gsl_meth = gsl_interp_polynomial; break;
        case OptionConstants::PRESCRIBED_DATA_INTERP_GSL_CSPLINE:
            gsl_meth = gsl_interp_cspline; break;
        case OptionConstants::PRESCRIBED_DATA_INTERP_GSL_AKIMA:
            gsl_meth = gsl_interp_akima; break;
        case OptionConstants::PRESCRIBED_DATA_INTERP_GSL_STEFFEN:
            gsl_meth = gsl_interp_steffen; break;

        default:
            throw SettingsException(
                "%s: Unrecognized interpolation method on radial grid: %d.",
                modname.c_str(), rinterp
            );
    }

    const len_t Nr_targ = rgrid->GetNr();
    gsl_interp *interp = gsl_interp_alloc(gsl_meth, nr_inp);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    // Construct new 'x' and 't' vectors (since Interpolator1D assumes
    // ownership of the data, and 'Settings' doesn't renounces its,
    // we allocate separate data for the 'Interpolator1D' object)
    real_t *new_x = new real_t[nt*Nr_targ];
    real_t *new_t = new real_t[nt];
    for (len_t it = 0; it < nt; it++) {
        gsl_interp_init(interp, r, x+(it*nr_inp), nr_inp);
        new_t[it] = t[it];

        for (len_t ir = 0; ir < Nr_targ; ir++) {
            real_t xr = rgrid->GetR(ir);
            new_x[it*Nr_targ + ir] = gsl_interp_eval(interp, r, x+(it*nr_inp), xr, acc);
        }

        gsl_interp_accel_reset(acc);
    }

    gsl_interp_accel_free(acc);
    gsl_interp_free(interp);
    
    return new FVM::Interpolator1D(nt, Nr_targ, new_t, new_x, interp1_meth);
}

