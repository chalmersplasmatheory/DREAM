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
#include "DREAM/MultiInterpolator1D.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"


using namespace DREAM;
using namespace std;


/**
 * Define options for a "radius+ion" 'data' section in
 * the specified module.
 */
void SimulationGenerator::DefineDataIonR(
    const string& modname, Settings *s,
    const string& name
) {
    const len_t ndim[2] = {0};

    s->DefineSetting(modname + "/" + name + "/r", "Radial grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/rinterp", "Interpolation method to use for radial grid interpolation.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_GSL_LINEAR);
    s->DefineSetting(modname + "/" + name + "/x", "Prescribed data.", 2, ndim, (real_t*)nullptr);
}

/**
 * Load a set of ion charge state radial density profiles
 * from the specified module and section of the settings.
 * The densities will be interpolated onto the given radial
 * grid.
 *
 * modname: Name of settings module from which to load the data.
 * rgrid:   Radial grid onto which to interpolate the profiles.
 * s:       Settings object to get data from.
 * nZ0:     Number of charge states expected.
 * name:    Name of variable containing the data.
 */
real_t *SimulationGenerator::LoadDataIonR(
    const string& modname, FVM::RadialGrid *rgrid, Settings *s,
    const len_t nZ0, const string& name
) {
    len_t xdims[2], nr_inp;

    const real_t *r = s->GetRealArray(modname + "/" + name + "/r", 1, &nr_inp);
    const real_t *x = s->GetRealArray(modname + "/" + name + "/x", 2, xdims);

    if (nZ0 != xdims[0] || nr_inp != xdims[1])
        throw SettingsException(
            "%s: Inconsistent dimensions of data. Data has "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " but "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " was expected.",
            (modname+"/"+name).c_str(), xdims[0], xdims[1], nZ0, nr_inp
        );

    // Check if dataset is empty
    if (nZ0 == 0 || nr_inp == 0)
        return nullptr;

    enum OptionConstants::prescribed_data_interp_gsl rinterp =
        (enum OptionConstants::prescribed_data_interp_gsl)s->GetInteger(modname + "/" + name + "/rinterp");

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

    return InterpolateIonR(
        rgrid, nr_inp, nZ0, r, x, gsl_meth
    );
}

/**
 * Interpolate a set of ion radial profiles from the given grid
 * to the grid defined by 'rgrid'.
 *
 * rgrid:    Radial grid to interpolate TO.
 * nr_inp:   Number of grid points in 'r' and 'x'.
 * nZ0:      Number of ion charge states.
 * r:        Radial grid on which the input data is defined.
 * x:        Data to interpolate.
 * gsl_meth: Interpolation method to use.
 */
real_t *SimulationGenerator::InterpolateIonR(
    FVM::RadialGrid *rgrid, const len_t nr_inp, const len_t nZ0,
    const real_t *r, const real_t *x, const gsl_interp_type *gsl_meth
) {
    const len_t Nr_targ = rgrid->GetNr();
    real_t *new_x = new real_t[nZ0*Nr_targ];

    // If the profile contains only one radial point, we assume
    // a uniform density profile
    if (nr_inp == 1) {
        for (len_t iZ = 0; iZ < nZ0; iZ++)
            for (len_t ir = 0; ir < Nr_targ; ir++)
                new_x[iZ*Nr_targ + ir] = x[iZ];
    } else {
        gsl_interp *interp = gsl_interp_alloc(gsl_meth, nr_inp);
        gsl_interp_accel *acc = gsl_interp_accel_alloc();

        // Construct a new 'x' vector and interpolate from the
        // input grid to the DREAM radial grid.
        for (len_t iZ = 0; iZ < nZ0; iZ++) {
            gsl_interp_init(interp, r, x+(iZ*nr_inp), nr_inp);

            for (len_t ir = 0; ir < Nr_targ; ir++) {
                real_t xr = rgrid->GetR(ir);

                // Out-of-bounds?
                if (xr < r[0]) {
                    // Extrapolate linearly!
                    real_t x0 = x[iZ*nr_inp+0], x1 = x[iZ*nr_inp+1];
                    real_t r0 = r[0], r1 = r[1];
                    real_t v  = x0 - (x1-x0)/(r1-r0)*(r0-xr);

                    new_x[iZ*Nr_targ + ir] = v > 0 ? v : 0;
                } else if (xr > r[nr_inp-1]) {
                    // Extrapolate linearly!
                    real_t x0 = x[iZ*nr_inp+nr_inp-2], x1 = x[iZ*nr_inp+nr_inp-1];
                    real_t r0 = r[nr_inp-2], r1 = r[nr_inp-1];
                    real_t v  = x1 + (x1-x0)/(r1-r0)*(xr-r1);

                    new_x[iZ*Nr_targ + ir] = v > 0 ? v : 0;
                } else
                    new_x[iZ*Nr_targ + ir] = gsl_interp_eval(interp, r, x+(iZ*nr_inp), xr, acc);
            }

            gsl_interp_accel_reset(acc);
        }

        gsl_interp_accel_free(acc);
        gsl_interp_free(interp);
    }

    return new_x;
}

/**
 * Define options for a "time+ion" 'data' section in the specified module.
 */
void SimulationGenerator::DefineDataIonT(
	const string& modname, Settings *s,
	const string& name
) {
	const len_t ndim[2] = {0};

	s->DefineSetting(modname + "/" + name + "/t", "Time grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
	s->DefineSetting(modname + "/" + name + "/tinterp", "Interpolation method to se for the time grid interpolation.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_LINEAR);
	s->DefineSetting(modname + "/" + name + "/x", "Prescribed data.", 2, ndim, (real_t*)nullptr);
}

/**
 * Load data from the 'data' section of the specified module. The data is
 * expected to depend on time and ion charge state.
 *
 * modname: Name of settings module from which to load the data.
 * s:       Settings object to get data from.
 * nZ0:     Number of charge states expected.
 * name:    Name of variable containing the data.
 */
MultiInterpolator1D *SimulationGenerator::LoadDataIonT(
	const string& modname, Settings *s, const len_t nZ0,
	const string& name
) {
	len_t xdims[2], nt;
    const real_t *t = s->GetRealArray(modname + "/" + name + "/t", 1, &nt);
    const real_t *x = s->GetRealArray(modname + "/" + name + "/x", 2, xdims);

	if (xdims[0] == 0 || xdims[1] == 0)
		return nullptr;

    if (nZ0 != xdims[0] || nt != xdims[1])
        throw SettingsException(
            "%s: Inconsistent dimensions of data. Data has "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " but "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " was expected.",
            (modname+"/"+name).c_str(), xdims[0], xdims[1], nZ0, nt
        );

    // Check if dataset is empty
    if (nZ0 == 0 || nt == 0)
        return nullptr;

    enum OptionConstants::prescribed_data_interp tinterp =
        (enum OptionConstants::prescribed_data_interp)s->GetInteger(modname + "/" + name + "/tinterp");

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

    // Construct new 'x' and 't' vectors (since Interpolator1D assumes
    // ownership of the data, and 'Settings' doesn't renounce its,
    // we allocate separate data for the 'Interpolator1D' object)
    real_t *new_x = new real_t[nZ0*nt];
    real_t *new_t = new real_t[nt];

    for (len_t it = 0; it < nt; it++)
        new_t[it] = t[it];

	for (len_t iZ = 0; iZ < nZ0; iZ++)
		for (len_t it = 0; it < nt; it++)
			new_x[iZ*nt + it] = x[iZ*nt + it];

    return new MultiInterpolator1D(
        nZ0, nt, 1, new_t, new_x, interp1_meth
    );
}

/**
 * Define options for a "radius+time+ion" 'data' section in
 * the specified module.
 */
void SimulationGenerator::DefineDataIonRT(
    const string& modname, Settings *s,
    const string& name
) {
    const len_t ndim[3] = {0};

    s->DefineSetting(modname + "/" + name + "/r", "Radial grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/rinterp", "Interpolation method to use for radial grid interpolation.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_GSL_LINEAR);
    s->DefineSetting(modname + "/" + name + "/t", "Time grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/tinterp", "Interpolation method to use for time grid interpolation.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_LINEAR);
    s->DefineSetting(modname + "/" + name + "/x", "Prescribed data.", 3, ndim, (real_t*)nullptr);
}

/**
 * Load data from the 'data' section of the specified module.
 * The data is expected to depend on radius, time and ion charge
 * state. It is interpolated in radius to the given radial grid.
 *
 * modname: Name of settings module from which to load the data.
 * rgrid:   Radial grid onto which to interpolate the profiles.
 * s:       Settings object to get data from.
 * nZ0:     Number of charge states expected.
 * name:    Name of variable containing the data.
 */
MultiInterpolator1D *SimulationGenerator::LoadDataIonRT(
    const string& modname, FVM::RadialGrid *rgrid, Settings *s,
    const len_t nZ0, const string& name
) {
    len_t xdims[3], nr_inp, nt;

    const real_t *r = s->GetRealArray(modname + "/" + name + "/r", 1, &nr_inp);
    const real_t *t = s->GetRealArray(modname + "/" + name + "/t", 1, &nt);
    const real_t *x = s->GetRealArray(modname + "/" + name + "/x", 3, xdims);

    if (nZ0 != xdims[0] || nt != xdims[1] || nr_inp != xdims[2])
        throw SettingsException(
            "%s: Inconsistent dimensions of data. Data has "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " but "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " was expected.",
            (modname+"/"+name).c_str(), xdims[0], xdims[1], xdims[2], nZ0, nt, nr_inp
        );

    // Check if dataset is empty
    if (nZ0 == 0 || nt == 0 || nr_inp == 0)
        return nullptr;

    enum OptionConstants::prescribed_data_interp_gsl rinterp =
        (enum OptionConstants::prescribed_data_interp_gsl)s->GetInteger(modname + "/" + name + "/rinterp");
    enum OptionConstants::prescribed_data_interp tinterp =
        (enum OptionConstants::prescribed_data_interp)s->GetInteger(modname + "/" + name + "/tinterp");

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

    // Construct new 'x' and 't' vectors (since Interpolator1D assumes
    // ownership of the data, and 'Settings' doesn't renounces its,
    // we allocate separate data for the 'Interpolator1D' object)
    real_t *new_x = new real_t[nZ0*nt*Nr_targ];
    real_t *new_t = new real_t[nt];

    for (len_t it = 0; it < nt; it++)
        new_t[it] = t[it];

    if (nr_inp == 1) {
        for (len_t iZ = 0; iZ < nZ0; iZ++)
            for (len_t it = 0; it < nt; it++)
                for (len_t ir = 0; ir < Nr_targ; ir++)
                    new_x[(iZ*nt + it)*Nr_targ + ir] = x[(iZ*nt + it)*nr_inp];
    } else {
        gsl_interp *interp = gsl_interp_alloc(gsl_meth, nr_inp);
        gsl_interp_accel *acc = gsl_interp_accel_alloc();

        for (len_t iZ = 0; iZ < nZ0; iZ++) {
            for (len_t it = 0; it < nt; it++) {
                gsl_interp_init(interp, r, x+((iZ*nt + it)*nr_inp), nr_inp);

                for (len_t ir = 0; ir < Nr_targ; ir++) {
                    real_t xr = rgrid->GetR(ir);

                    if (xr < r[0]) {
                        // Extrapolate linearly!
                        real_t x0 = x[(iZ*nt + it)*nr_inp+0], x1 = x[(iZ*nt + it)*nr_inp+1];
                        real_t r0 = r[0], r1 = r[1];
                        real_t v  = x0 - (x1-x0)/(r1-r0)*(r0-xr);

                        new_x[(iZ*nt + it)*Nr_targ + ir] = v > 0 ? v : 0;
                    } else if (xr > r[nr_inp-1]) {
                        // Extrapolate linearly!
                        real_t x0 = x[(iZ*nt + it)*nr_inp+nr_inp-2], x1 = x[(iZ*nt + it)*nr_inp+nr_inp-1];
                        real_t r0 = r[nr_inp-2], r1 = r[nr_inp-1];
                        real_t v  = x1 + (x1-x0)/(r1-r0)*(xr-r1);

                        new_x[(iZ*nt + it)*Nr_targ + ir] = v > 0 ? v : 0;
                    } else
                        new_x[(iZ*nt + it)*Nr_targ + ir] = gsl_interp_eval(interp, r, x+((iZ*nt + it)*nr_inp), xr, acc);
                }

                gsl_interp_accel_reset(acc);
            }
        }

        gsl_interp_accel_free(acc);
        gsl_interp_free(interp);
    }

    return new MultiInterpolator1D(
        nZ0, nt, Nr_targ, new_t, new_x, interp1_meth
    );
}

/**
 * Define options for a radial profile data section in the
 * specified module.
 */
void SimulationGenerator::DefineDataR(
    const string& modname, Settings *s,
    const string& name
) {
    s->DefineSetting(modname + "/" + name + "/rinterp", "Interpolation method to use for radial grid interpolation.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_GSL_LINEAR);
    s->DefineSetting(modname + "/" + name + "/r", "Radial grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/x", "Prescribed data.", 0, (real_t*)nullptr);
}

/**
 * Load data from the 'data' section of the specified module.
 * The data is expected to depend only on radius. It is
 * interpolated in radius to the given radial grid.
 *
 * modname: Name of the module to load data from.
 * rgrid:   Radial grid to inerpolate data to.
 * s:       Settings object to load data from.
 * name:    Name of group containing the data structure (default: "data").
 */
real_t *SimulationGenerator::LoadDataR(
    const string& modname, FVM::RadialGrid *rgrid, Settings *s,
    const string& name
) {
    len_t nx, nr_inp;

    const real_t *r = s->GetRealArray(modname + "/" + name + "/r", 1, &nr_inp);
    const real_t *x = s->GetRealArray(modname + "/" + name + "/x", 1, &nx);

    if (nr_inp != nx)
        throw SettingsException(
            "%s: Inconsistent size of data. Data has "
            LEN_T_PRINTF_FMT " elements, but " LEN_T_PRINTF_FMT " were expected.",
            (modname+"/"+name).c_str(), nx, nr_inp
        );

    enum OptionConstants::prescribed_data_interp_gsl rinterp =
        (enum OptionConstants::prescribed_data_interp_gsl)s->GetInteger(modname + "/" + name + "/rinterp");

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

    // Interpolate given profile to computational grid
    if (nr_inp == 0)
        return nullptr;
    else
        return SimulationGenerator::InterpolateR(nr_inp, r, x, rgrid, gsl_meth);
}

/**
 * Interpolate a radial profile from the given grid to the
 * grid defined by 'rgrid'.
 *
 * nr_inp:   Number of grid points in 'r' and 'x'.
 * r:        Radial grid on which the input data is defined.
 * x:        Data to interpolate.
 * rgrid:    Radial grid to interpolate TO.
 * gsl_meth: Interpolation method to use.
 */
real_t *SimulationGenerator::InterpolateR(
    const len_t nr_inp, const real_t *r, const real_t *x,
    FVM::RadialGrid *rgrid, const gsl_interp_type *gsl_meth
) {
    const len_t Nr_targ = rgrid->GetNr();
    real_t *new_x = new real_t[Nr_targ];

    if (nr_inp == 1) {
        // Uniform profile
        for (len_t ir = 0; ir < Nr_targ; ir++)
            new_x[ir] = x[0];
    } else {
        gsl_interp *interp = gsl_interp_alloc(gsl_meth, nr_inp);
        gsl_interp_accel *acc = gsl_interp_accel_alloc();

        gsl_interp_init(interp, r, x, nr_inp);

        for (len_t ir = 0; ir < Nr_targ; ir++) {
            real_t xr = rgrid->GetR(ir);
            // Out-of-bounds?
            if (xr < r[0]) {
                // Extrapolate linearly!
                real_t x0 = x[0], x1 = x[1];
                real_t r0 = r[0], r1 = r[1];
                real_t v  = x0 - (x1-x0)/(r1-r0)*(r0-xr);

                new_x[ir] = v > 0 ? v : 0;
            } else if (xr > r[nr_inp-1]) {
                // Extrapolate linearly!
                real_t x0 = x[nr_inp-2], x1 = x[nr_inp-1];
                real_t r0 = r[nr_inp-2], r1 = r[nr_inp-1];
                real_t v  = x1 + (x1-x0)/(r1-r0)*(xr-r1);

                new_x[ir] = v > 0 ? v : 0;
            } else
                new_x[ir] = gsl_interp_eval(interp, r, x, xr, acc);
        }

        gsl_interp_accel_free(acc);
        gsl_interp_free(interp);
    }

    return new_x;
}

/**
 * Define options for a temporal 'data' section in the
 * specified module.
 */
void SimulationGenerator::DefineDataT(
    const string& modname, Settings *s, const string& name
) {
    s->DefineSetting(modname + "/" + name + "/t", "Time grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/tinterp", "Interpolation method to use for time grid interpolation.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_LINEAR);
    s->DefineSetting(modname + "/" + name + "/x", "Prescribed data.", 0, (real_t*)nullptr);
}

/**
 * Load data from the 'data' seciton of the specified module.
 * The data is expected to depend on time only.
 *
 * modname: Name of module to load data from.
 * s:       Settings object to load data from.
 * name:    Name of group containing data structure (default: "data").
 */
FVM::Interpolator1D *SimulationGenerator::LoadDataT(
    const string& modname, Settings *s, const string& name
) {
    len_t nx, nt;

    const real_t *t = s->GetRealArray(modname + "/" + name + "/t", 1, &nt);
    const real_t *x = s->GetRealArray(modname + "/" + name + "/x", 1, &nx);

    if (nt != nx)
        throw SettingsException(
            "%s: Inconsistent dimensions of data. Data has "
            LEN_T_PRINTF_FMT " elements while the time vector has "
            LEN_T_PRINTF_FMT " elements.",
            (modname+"/"+name).c_str(), nx, nt
        );
	
	if (nt == 0)
		return nullptr;

    enum OptionConstants::prescribed_data_interp tinterp =
        (enum OptionConstants::prescribed_data_interp)s->GetInteger(modname + "/" + name + "/tinterp");

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

    real_t *new_x = new real_t[nt];
    real_t *new_t = new real_t[nt];

    // Copy data
    for (len_t it = 0; it < nt; it++) {
        new_t[it] = t[it];
        new_x[it] = x[it];
    }
    
    return new FVM::Interpolator1D(nt, 1, new_t, new_x, interp1_meth);
}

/**
 * Define options for a "radius+time" 'data' section in the
 * specified module.
 */
void SimulationGenerator::DefineDataRT(
    const string& modname, Settings *s,
    const string& name
) {
    const len_t ndim[2] = {0,0};

    s->DefineSetting(modname + "/" + name + "/r", "Radial grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/rinterp", "Interpolation method to use for radial grid interpolation.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_GSL_LINEAR);
    s->DefineSetting(modname + "/" + name + "/t", "Time grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/tinterp", "Interpolation method to use for time grid interpolation.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_LINEAR);
    s->DefineSetting(modname + "/" + name + "/x", "Prescribed data.", 2, ndim, (real_t*)nullptr);
}

/**
 * Load data from the 'data' section of the specified module.
 * The data is expected to depend on radius and time. It is
 * interpolated in radius to the given radial grid.
 *
 * modname:   Name of module to load data from.
 * rgrid:     Radial grid to interpolate data to.
 * s:         Settings object to load data from.
 * name:      Name of group containing data structure (default: "data").
 * rFluxGrid: If 'true', interpolates the loaded quantity to the
 *            radial flux grid, instead of the distribution grid.
 */
struct dream_2d_data *SimulationGenerator::LoadDataRT(
    const string& modname, FVM::RadialGrid *rgrid, Settings *s,
    const string& name, const bool rFluxGrid
) {
    len_t xdims[2], nr_inp, nt;

    const real_t *r = s->GetRealArray(modname + "/" + name + "/r", 1, &nr_inp);
    const real_t *t = s->GetRealArray(modname + "/" + name + "/t", 1, &nt);
    const real_t *x = s->GetRealArray(modname + "/" + name + "/x", 2, xdims);

    if (nt != xdims[0] || nr_inp != xdims[1])
        throw SettingsException(
            "%s: Inconsistent dimensions of data. Data has "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " but "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " was expected.",
            (modname+"/"+name).c_str(), xdims[0], xdims[1], nt, nr_inp
        );

    enum OptionConstants::prescribed_data_interp_gsl rinterp =
        (enum OptionConstants::prescribed_data_interp_gsl)s->GetInteger(modname + "/" + name + "/rinterp");
    enum OptionConstants::prescribed_data_interp tinterp =
        (enum OptionConstants::prescribed_data_interp)s->GetInteger(modname + "/" + name + "/tinterp");

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

    const len_t Nr_targ = rgrid->GetNr() + rFluxGrid;

    real_t *new_x = new real_t[nt*Nr_targ];
    real_t *new_t = new real_t[nt];

    if (nr_inp == 1) {
        // Uniform profile
        for (len_t it = 0; it < nt; it++) {
            new_t[it] = t[it];
            for (len_t ir = 0; ir < Nr_targ; ir++)
                new_x[it*Nr_targ + ir] = x[it];
        }
    } else {
        gsl_interp *interp = gsl_interp_alloc(gsl_meth, nr_inp);
        gsl_interp_accel *acc = gsl_interp_accel_alloc();

        // Construct new 'x' and 't' vectors (since Interpolator1D assumes
        // ownership of the data, and 'Settings' doesn't renounces its,
        // we allocate separate data for the 'Interpolator1D' object)
        for (len_t it = 0; it < nt; it++) {
            gsl_interp_init(interp, r, x+(it*nr_inp), nr_inp);
            new_t[it] = t[it];

            for (len_t ir = 0; ir < Nr_targ; ir++) {
                real_t xr;
                if (rFluxGrid)
                    xr = rgrid->GetR_f(ir);
                else
                    xr = rgrid->GetR(ir);

                if (xr < r[0]) {
                    // Extrapolate linearly!
                    real_t x0 = x[it*nr_inp + 0], x1 = x[it*nr_inp + 1];
                    real_t r0 = r[0], r1 = r[1];
                    real_t v  = x0 - (x1-x0)/(r1-r0)*(r0-xr);

                    new_x[it*Nr_targ + ir] = v > 0 ? v : 0;
                } else if (xr > r[nr_inp-1]) {
                    // Extrapolate linearly!
                    real_t x0 = x[it*nr_inp+nr_inp-2], x1 = x[it*nr_inp+nr_inp-1];
                    real_t r0 = r[nr_inp-2], r1 = r[nr_inp-1];
                    real_t v  = x1 + (x1-x0)/(r1-r0)*(xr-r1);

                    new_x[it*Nr_targ + ir] = v > 0 ? v : 0;
                } else
                    new_x[it*Nr_targ + ir] = gsl_interp_eval(interp, r, x+(it*nr_inp), xr, acc);
            }

            gsl_interp_accel_reset(acc);
        }

        gsl_interp_accel_free(acc);
        gsl_interp_free(interp);
    }

    struct dream_2d_data *d = new struct dream_2d_data;

    d->x = new_x;
    d->t = new_t;
    d->r = new real_t[Nr_targ];
    for (len_t i = 0; i < Nr_targ; i++)
        if (rFluxGrid)
            d->r[i] = rgrid->GetR_f(i);
        else
            d->r[i] = rgrid->GetR(i);
    d->nt = nt;
    d->nr = Nr_targ;
    d->interp = interp1_meth;

    return d;
}
    
/**
 * Load data from the 'data' section of the specified module.
 * The data is expected to depend on radius and time. It is
 * interpolated in radius to the given radial grid.
 *
 * modname:   Name of module to load data from.
 * rgrid:     Radial grid to interpolate data to.
 * s:         Settings object to load data from.
 * name:      Name of group containing data structure (default: "data").
 * rFluxGrid: If 'true', interpolates the loaded quantity to the
 *            radial flux grid, instead of the distribution grid.
 */
FVM::Interpolator1D *SimulationGenerator::LoadDataRT_intp(
    const string& modname, FVM::RadialGrid *rgrid, Settings *s,
    const string& name, const bool rFluxGrid
) {
    struct dream_2d_data *d = LoadDataRT(modname, rgrid, s, name, rFluxGrid);

    auto i = new FVM::Interpolator1D(d->nt, d->nr, d->t, d->x, d->interp);

    // d->r is the only array not used here anymore...
    delete [] d->r;
    // Also delete 'd', since 'i' only needs the actual data...
    delete d;

    return i;
}

/**
 * Define options for a "radius+momentum+momentum" 'data'
 * section in the specified module.
 */
void SimulationGenerator::DefineDataR2P(
    const string& modname, Settings *s,
    const string& name
) {
    const len_t ndim[3] = {0,0,0};

    s->DefineSetting(modname + "/" + name + "/interp", "3D interpolation method to use.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP3D_LINEAR);
    s->DefineSetting(modname + "/" + name + "/r", "Radial grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/p", "Momentum grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/ppar", "Parallel momentum grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/pperp", "Perpendicular momentum grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/x", "Prescribed data.", 3, ndim, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/xi", "Pitch grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
}

/**
 * Load data from the 'data' section of the specified module.
 * The data is expected to depend on radius and two momentum
 * coordinates. It is interpolated to the given grid.
 *
 * modname: Name of module to load data from.
 * grid:    Grid to interpolate data to.
 * s:       Settings object to load data from.
 * name:    Name of group containing data structure (default: "data").
 */
FVM::Interpolator3D *SimulationGenerator::LoadDataR2P(
    const string& modname, Settings *s, const string& name
) {
    len_t xdims[3], nr, np1, np2;

    const real_t *_r = s->GetRealArray(modname + "/" + name + "/r", 1, &nr);
    const real_t *_x = s->GetRealArray(modname + "/" + name + "/x", 3, xdims);

    enum OptionConstants::prescribed_data_interp3d meth =
        (enum OptionConstants::prescribed_data_interp3d)s->GetInteger(modname + "/" + name + "/interp");

    // Select Interpolator3D interpolation method
    enum FVM::Interpolator3D::interp_method interp_meth;
    switch (meth) {
        case OptionConstants::PRESCRIBED_DATA_INTERP3D_NEAREST:
            interp_meth = FVM::Interpolator3D::INTERP_NEAREST; break;
        case OptionConstants::PRESCRIBED_DATA_INTERP3D_LINEAR:
            interp_meth = FVM::Interpolator3D::INTERP_LINEAR; break;

        default:
            throw SettingsException(
                "%s: Unrecognized 3D interpolation method: %d.",
                modname.c_str(), meth
            );
    }

    // Load momentum grid vectors
    const real_t *_p1, *_p2;
    FVM::Interpolator3D::momentumgrid_type momtype;

    if ((_p1=s->GetRealArray(modname + "/" + name + "/p", 1, &np1, false)) != nullptr &&
        (_p2=s->GetRealArray(modname + "/" + name + "/xi", 1, &np2, false)) != nullptr) {

        momtype = FVM::Interpolator3D::GRID_PXI;
		s->MarkUsed(modname + "/" + name + "/p");
		s->MarkUsed(modname + "/" + name + "/xi");
    } else if ((_p1=s->GetRealArray(modname + "/" + name + "/ppar", 1, &np1, false)) != nullptr &&
        (_p2=s->GetRealArray(modname + "/" + name + "/pperp", 1, &np2, false)) != nullptr) {

        momtype = FVM::Interpolator3D::GRID_PPARPPERP;
		s->MarkUsed(modname + "/" + name + "/ppar");
		s->MarkUsed(modname + "/" + name + "/pperp");
    } else
        throw SettingsException(
            "%s: No momentum grid set for data.",
            modname.c_str()
        );

    // Verify array lengths and dimensions
    if (xdims[0] != nr || xdims[1] != np2 || xdims[2] != np1)
        throw SettingsException(
            "%s: Invalid dimensions of data: " LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT
            ". Expected: " LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT ".",
            xdims[0], xdims[1], xdims[2], nr, np2, np1
        );

    // Copy data
    real_t *x  = new real_t[nr*np1*np2], *r = new real_t[nr];
    real_t *p1 = new real_t[np1];
    real_t *p2 = new real_t[np2];

    for (len_t i = 0; i < nr*np1*np2; i++)
        x[i] = _x[i];
    for (len_t i = 0; i < nr; i++)
        r[i] = _r[i];
    for (len_t i = 0; i < np1; i++)
        p1[i] = _p1[i];
    for (len_t i = 0; i < np2; i++)
        p2[i] = _p2[i];

    // Finally, construct Interpolator3D object
    FVM::Interpolator3D *interp = new FVM::Interpolator3D(
        nr, np2, np1, r, p2, p1, x,
        momtype, interp_meth
    );

    return interp;
}

/**
 * Define options for a "time+radius+momentum+momentum" 'data'
 * section in the specified module.
 */
void SimulationGenerator::DefineDataTR2P(
    const string& modname, Settings *s,
    const string& name
) {
    const len_t ndim[4] = {0,0,0,0};

    s->DefineSetting(modname + "/" + name + "/interp3d", "3D interpolation method to use.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP3D_LINEAR);
    s->DefineSetting(modname + "/" + name + "/interp1d", "Time interpolation method to use.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_LINEAR);
    s->DefineSetting(modname + "/" + name + "/p", "Momentum grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/ppar", "Parallel momentum grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/pperp", "Perpendicular momentum grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/r", "Radial grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/t", "Time grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/x", "Prescribed data.", 4, ndim, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/xi", "Pitch grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
}

/**
 * Load data from the 'data' section of the specified module.
 * The data is expected to depend on radius and two momentum
 * coordinates. It is interpolated to the given grid.
 *
 * modname: Name of module to load data from.
 * s:       Settings object to load data from.
 * name:    Name of group containing data structure (default: "data").
 */
struct dream_4d_data *SimulationGenerator::LoadDataTR2P(
    const string& modname, Settings *s, const string& name
) {
    len_t xdims[4], nt, nr, np1, np2;

    const real_t *_r = s->GetRealArray(modname + "/" + name + "/r", 1, &nr);
    const real_t *_t = s->GetRealArray(modname + "/" + name + "/t", 1, &nt);
    const real_t *_x = s->GetRealArray(modname + "/" + name + "/x", 4, xdims);

    enum OptionConstants::prescribed_data_interp3d meth3d =
        (enum OptionConstants::prescribed_data_interp3d)s->GetInteger(modname + "/" + name + "/interp3d");
    enum OptionConstants::prescribed_data_interp meth1d =
        (enum OptionConstants::prescribed_data_interp)s->GetInteger(modname + "/" + name + "/interp1d");

    // Select Interpolator3D interpolation method
    enum FVM::Interpolator3D::interp_method interp3d;
    switch (meth3d) {
        case OptionConstants::PRESCRIBED_DATA_INTERP3D_NEAREST:
            interp3d = FVM::Interpolator3D::INTERP_NEAREST; break;
        case OptionConstants::PRESCRIBED_DATA_INTERP3D_LINEAR:
            interp3d = FVM::Interpolator3D::INTERP_LINEAR; break;

        default:
            throw SettingsException(
                "%s: Unrecognized 3D interpolation method: %d.",
                modname.c_str(), interp3d
            );
    }

    // Select Interpolator1D interpolation method
    enum FVM::Interpolator1D::interp_method interp1d;
    switch (meth1d) {
        case OptionConstants::PRESCRIBED_DATA_INTERP_NEAREST:
            interp1d = FVM::Interpolator1D::INTERP_NEAREST; break;
        case OptionConstants::PRESCRIBED_DATA_INTERP_LINEAR:
            interp1d = FVM::Interpolator1D::INTERP_LINEAR; break;

        default:
            throw SettingsException(
                "%s: Unrecognized 1D interpolation method: %d.",
                modname.c_str(), interp1d
            );
    }

    // Load momentum grid vectors
    const real_t *_p1, *_p2;
    FVM::Interpolator3D::momentumgrid_type momtype;

    if ((_p1=s->GetRealArray(modname + "/" + name + "/p", 1, &np1, false)) != nullptr &&
        (_p2=s->GetRealArray(modname + "/" + name + "/xi", 1, &np2, false)) != nullptr) {

        momtype = FVM::Interpolator3D::GRID_PXI;
		s->MarkUsed(modname + "/" + name + "/p");
		s->MarkUsed(modname + "/" + name + "/xi");
    } else if ((_p1=s->GetRealArray(modname + "/" + name + "/ppar", 1, &np1, false)) != nullptr &&
        (_p2=s->GetRealArray(modname + "/" + name + "/pperp", 1, &np2, false)) != nullptr) {

        momtype = FVM::Interpolator3D::GRID_PPARPPERP;
		s->MarkUsed(modname + "/" + name + "/ppar");
		s->MarkUsed(modname + "/" + name + "/pperp");
    } else
        throw SettingsException(
            "%s: No momentum grid set for data.",
            modname.c_str()
        );

    // Verify dimensions...
    if (xdims[0] != nt || xdims[1] != nr || xdims[2] != np2 || xdims[3] != np1)
        throw SettingsException(
            "%s: Invalid dimensions of data: " LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT
            ". Expected: " LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT ".",
            xdims[0], xdims[1], xdims[2], xdims[3], nt, nr, np2, np1
        );

    // Copy data...
    real_t **x = new real_t*[nt];
    x[0] = new real_t[nt*nr*np1*np2];
    real_t *t = new real_t[nt];
    real_t *r = new real_t[nr];
    real_t *p1 = new real_t[np1];
    real_t *p2 = new real_t[np2];

    const len_t N = nr*np1*np2;
    for (len_t j = 0; j < nt; j++) {
        if (j > 0)
            x[j] = x[j-1] + N;

        for (len_t i = 0; i < N; i++)
            x[j][i] = _x[i + j*N];
    }
    for (len_t i = 0; i < nt; i++)
        t[i] = _t[i];
    for (len_t i = 0; i < nr; i++)
        r[i] = _r[i];
    for (len_t i = 0; i < np1; i++)
        p1[i] = _p1[i];
    for (len_t i = 0; i < np2; i++)
        p2[i] = _p2[i];


    // Construct the struct to return...
    struct dream_4d_data *d = new struct dream_4d_data;

    d->x   = x;
    d->t   = t;
    d->r   = r;
    d->p1  = p1;
    d->p2  = p2;
    d->nt  = nt;
    d->nr  = nr;
    d->np1 = np1;
    d->np2 = np2;

    d->gridtype = momtype;
    d->ps_interp = interp3d;
    d->time_interp = interp1d;

    return d;
}

