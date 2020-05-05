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
#include "DREAM/IonInterpolator1D.hpp"
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

    real_t *r = s->GetRealArray(modname + "/" + name + "/r", 1, &nr_inp);
    real_t *x = s->GetRealArray(modname + "/" + name + "/x", 2, xdims);

    if (nZ0 != xdims[0] || nr_inp != xdims[1])
        throw SettingsException(
            "%s: Inconsistent dimensions of data. Data has "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " but "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " was expected.",
            xdims[0], xdims[1], nZ0, nr_inp
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

    const len_t Nr_targ = rgrid->GetNr();
    gsl_interp *interp = gsl_interp_alloc(gsl_meth, nr_inp);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    // Construct a new 'x' vector and interpolate from the
    // input grid to the DREAM radial grid.
    real_t *new_x = new real_t[nZ0*Nr_targ];
    for (len_t iZ = 0; iZ < nZ0; iZ++) {
        gsl_interp_init(interp, r, x+(iZ*nr_inp), nr_inp);

        for (len_t ir = 0; ir < Nr_targ; ir++) {
            real_t xr = rgrid->GetR(ir);
            new_x[iZ*Nr_targ + ir] = gsl_interp_eval(interp, r, x+(iZ*nr_inp), xr, acc);
        }

        gsl_interp_accel_reset(acc);
    }

    gsl_interp_accel_free(acc);
    gsl_interp_free(interp);

    return new_x;
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
    s->DefineSetting(modname + "/" + name + "/tinterp", "Interpolation method to use for time grid interpolation.", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_GSL_LINEAR);
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
IonInterpolator1D *SimulationGenerator::LoadDataIonRT(
    const string& modname, FVM::RadialGrid *rgrid, Settings *s,
    const len_t nZ0, const string& name
) {
    len_t xdims[3], nr_inp, nt;

    real_t *r = s->GetRealArray(modname + "/" + name + "/r", 1, &nr_inp);
    real_t *t = s->GetRealArray(modname + "/" + name + "/t", 1, &nt);
    real_t *x = s->GetRealArray(modname + "/" + name + "/x", 3, xdims);

    if (nZ0 != xdims[0] || nt != xdims[1] || nr_inp != xdims[2])
        throw SettingsException(
            "%s: Inconsistent dimensions of data. Data has "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " but "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " was expected.",
            xdims[0], xdims[1], xdims[2], nZ0, nt, nr_inp
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

    const len_t Nr_targ = rgrid->GetNr();
    gsl_interp *interp = gsl_interp_alloc(gsl_meth, nr_inp);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    // Construct new 'x' and 't' vectors (since Interpolator1D assumes
    // ownership of the data, and 'Settings' doesn't renounces its,
    // we allocate separate data for the 'Interpolator1D' object)
    real_t *new_x = new real_t[nZ0*nt*Nr_targ];
    real_t *new_t = new real_t[nt];

    for (len_t it = 0; it < nt; it++)
        new_t[it] = t[it];

    for (len_t iZ = 0; iZ < nZ0; iZ++) {
        for (len_t it = 0; it < nt; it++) {
            gsl_interp_init(interp, r, x+((iZ*nt + it)*nr_inp), nr_inp);

            for (len_t ir = 0; ir < Nr_targ; ir++) {
                real_t xr = rgrid->GetR(ir);
                new_x[(iZ*nt + it)*Nr_targ + ir] = gsl_interp_eval(interp, r, x+((iZ*nt + it)*nr_inp), xr, acc);
            }

            gsl_interp_accel_reset(acc);
        }
    }

    gsl_interp_accel_free(acc);
    gsl_interp_free(interp);

    return new IonInterpolator1D(
        nZ0, nt, Nr_targ, new_t, new_x, interp1_meth
    );
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
 * modname: Name of module to load data from.
 * rgrid:   Radial grid to interpolate data to.
 * s:       Settings object to load data from.
 * name:    Name of group containing data structure (default: "data").
 */
FVM::Interpolator1D *SimulationGenerator::LoadDataRT(
    const string& modname, FVM::RadialGrid *rgrid, Settings *s,
    const string& name
) {
    len_t xdims[2], nr_inp, nt;

    real_t *r = s->GetRealArray(modname + "/" + name + "/r", 1, &nr_inp);
    real_t *t = s->GetRealArray(modname + "/" + name + "/t", 1, &nt);
    real_t *x = s->GetRealArray(modname + "/" + name + "/x", 2, xdims);

    if (nt != xdims[0] || nr_inp != xdims[1])
        throw SettingsException(
            "%s: Inconsistent dimensions of data. Data has "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " but "
            LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " was expected.",
            xdims[0], xdims[1], nt, nr_inp
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

    real_t *_r = s->GetRealArray(modname + "/" + name + "/r", 1, &nr);
    real_t *_x = s->GetRealArray(modname + "/" + name + "/x", 3, xdims);

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
    real_t *_p1, *_p2;
    FVM::Interpolator3D::momentumgrid_type momtype;

    if ((_p1=s->GetRealArray(modname + "/" + name + "/p", 1, &np1, false)) != nullptr &&
        (_p2=s->GetRealArray(modname + "/" + name + "/xi", 1, &np2, false)) != nullptr) {

        momtype = FVM::Interpolator3D::GRID_PXI;
    } else if ((_p1=s->GetRealArray(modname + "/" + name + "/ppar", 1, &np1, false)) != nullptr &&
        (_p2=s->GetRealArray(modname + "/" + name + "/pperp", 1, &np2, false)) != nullptr) {

        momtype = FVM::Interpolator3D::GRID_PPARPPERP;
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

