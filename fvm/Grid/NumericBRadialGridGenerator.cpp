/**
 * Implementation of the numeric magnetic field radial grid generator. This
 * grid generator loads a numeric magnetic field from the specified file and
 * builds a correspondingly shaped radial grid.
 */

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline2d.h>
#include "FVM/Grid/NumericBRadialGridGenerator.hpp"


using namespace DREAM::FVM;


/**
 * Constructor for a uniform radial grid.
 *
 * nr: Number of radial grid points in uniform radial *distribution* grid.
 * r0: Value of innermost point on radial *flux* grid.
 * ra: Value of outermost point on radial *flux* grid.
 * mf: Name of file containing the magnetic field data to load.
 */
NumericBRadialGridGenerator::NumericBRadialGridGenerator(
    const len_t nr, const real_t r0, const real_t ra,
    const std::string& mf
) : RadialGridGenerator(nr), rMin(r0), rMax(ra) {

    this->isUpDownSymmetric = false;
    LoadMagneticFieldData(mf);
}

/**
 * Constructor for a custom radial grid.
 *
 * r_f: Radial flux grid desired.
 * nr:  Number of points on distribution grids corresponding to 
 *      'r_f' (i.e. the number of points in 'r_f', minus one).
 */
NumericBRadialGridGenerator::NumericBRadialGridGenerator(
    const real_t *r_f, const len_t nr,
    const std::string& mf
) : RadialGridGenerator(nr) {

    this->isUpDownSymmetric = false;

    this->rf_provided = new real_t[nr];
    for (len_t i = 0; i < nr+1; i++)
        this->rf_provided[i] = r_f[i];

    LoadMagneticFieldData(mf);
}

/**
 * Destructor.
 */
NumericBRadialGridGenerator::~NumericBRadialGridGenerator() {
    if (this->rf_provided != nullptr)
        delete [] this->rf_provided;

    if (this->spline_R != nullptr) {
        gsl_spline2d_free(this->spline_R);
        gsl_spline2d_free(this->spline_Z);
        gsl_spline2d_free(this->spline_BR);
        gsl_spline2d_free(this->spline_BZ);
        gsl_spline2d_free(this->spline_Bphi);
    }
}


/**
 * Load magnetic field data from the named file.
 *
 * filename: Name of file to load data from.
 */
void NumericBRadialGridGenerator::LoadMagneticFieldData(
    const std::string& filename
) {
    SFile *sf = SFile::Create(filename, SFILE_MODE_READ);
    this->LoadMagneticFieldData(sf);

    sf->Close();
    delete sf;
}

/**
 * Load magnetic field data from the given file.
 *
 * sf: SFile object representing the file.
 */
void NumericBRadialGridGenerator::LoadMagneticFieldData(
    SFile *sf
) {
    sfilesize_t fsize[2];

    #define ASSERT_DIMS(var) \
        if (fsize[0] != this->ntheta || fsize[1] != this->npsi) \
            throw FVMException( \
                "%s: Invalid dimensions of vector '" var "' (%llu, %llu). " \
                "Expected (" LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT ").", \
                sf->filename.c_str(), fsize[0], fsize[1], ntheta, npsi \
            )

    this->name = sf->GetString("equil/id");

    // Magnetic axis coordinates
    this->Rp = (real_t)sf->GetScalar("equil/Rp");
    this->Zp = (real_t)sf->GetScalar("equil/Zp");

    // Poloidal flux coordinate grid
    double *_psi = sf->GetList("equil/psi_apRp", fsize);
    this->npsi = fsize[1]==1 ? fsize[0] : fsize[1];

    // Poloidal angle coordinate grid
    double *_theta = sf->GetList("equil/theta", fsize);
    this->ntheta = fsize[1]==1 ? fsize[0] : fsize[1];

    // Radial meshgrid
    double *_R = sf->GetList("equil/ptx", fsize); ASSERT_DIMS("ptx");
    double *_Z = sf->GetList("equil/pty", fsize); ASSERT_DIMS("pty");

    double *_Br = sf->GetList("equil/ptBx", fsize); ASSERT_DIMS("ptBx");
    double *_Bz = sf->GetList("equil/ptBy", fsize); ASSERT_DIMS("ptBz");
    double *_Bp = sf->GetList("equil/ptBPHI", fsize); ASSERT_DIMS("ptBPHI");

    auto convert_data = [this](const double *a, const len_t nx, const len_t ny) {
        if (typeid(real_t) == typeid(double))
            return (real_t*)a;

        real_t *d = new real_t[nx*ny];
        for (len_t i = 0; i < nx*ny; i++)
            d[i] = (real_t)a[i];

        delete [] a;
        return d;
    };

    // Set magnetic field data
    this->psi      = convert_data(_psi, npsi, 1);
    this->theta    = convert_data(_theta, ntheta, 1);
    
    this->R        = convert_data(_R, ntheta, npsi);
    this->Z        = convert_data(_Z, ntheta, npsi);

    this->dataBR   = convert_data(_Br, ntheta, npsi);
    this->dataBZ   = convert_data(_Bz, ntheta, npsi);
    this->dataBphi = convert_data(_Bp, ntheta, npsi);

    // Evaluate minor radius in outer midplane
    this->input_r = new real_t[this->npsi];
    for (len_t i = 0; i < this->npsi; i++)
        this->input_r[i] = R[i];
}

/**
 * Rebuild this numeric radial B grid.
 */
bool NumericBRadialGridGenerator::Rebuild(const real_t, RadialGrid *rGrid) {
    this->r   = new real_t[GetNr()];
    this->r_f = new real_t[GetNr()+1];

    real_t
        *dr   = new real_t[GetNr()],
        *dr_f = new real_t[GetNr()-1];

    // If a custom grid has been specified, set it here,
    // otherwise we generate a uniform radial grid
    if (rf_provided == nullptr) {
        if (rMin < 0)
            throw FVMException("NumericBRadialGrid: rMin < 0.");
        else if (rMax > this->input_r[this->npsi-1])
            throw FVMException("NumericBRadialGrid: Maximum r available in numeric magnetic field data is rMax = %.3f.", this->input_r[this->npsi-1]);
        else if (rMin >= rMax)
            throw FVMException("NumericBRadialGrid: rMin must be strictly less than rMax.");

        for (len_t i = 0; i < GetNr(); i++)
            dr[i] = (rMax-rMin) / GetNr();
        for (len_t i = 0; i < GetNr()+1; i++)
            r_f[i] = rMin + i*dr[0];
    } else {
        if (r_f[0] < 0)
            throw FVMException("NumericBRadialGrid: First point on custom radial grid is less than zero.");
        else if (r_f[GetNr()] > this->input_r[this->npsi-1])
            throw FVMException(
                "NumericBRadialGrid: Last point on custom radial grid may not be greater than "
                "the maximum r available in the numeric magnetic field data, rMax = %.3f.",
                this->input_r[this->npsi-1]
            );
        else if (r_f[0] >= r_f[GetNr()])
            throw FVMException("NumericBRadialGrid: The first point on the custom radial grid must be strictly less than the last point.");

        for (len_t i = 0; i < GetNr()+1; i++)
            r_f[i] = rf_provided[i];
        for (len_t i = 0; i < GetNr(); i++)
            dr[i] = r_f[i+1] - r_f[i];
    }

    // On the next rebuild, we would rather like to create
    // a new uniform grid at the newly set resolution
    delete [] rf_provided;
    rf_provided = nullptr;

    // Construct cell grid
    for (len_t i = 0; i < GetNr(); i++)
        r[i] = 0.5 * (r_f[i+1] + r_f[i]);
    for (len_t i = 0; i < GetNr()-1; i++)
        dr_f[i] = r[i+1] - r[i];

    rGrid->Initialize(r, r_f, dr, dr_f);

    // Construct splines for input data
    this->spline_R    = gsl_spline2d_alloc(gsl_interp2d_bilinear, this->npsi, this->ntheta);
    this->spline_Z    = gsl_spline2d_alloc(gsl_interp2d_bilinear, this->npsi, this->ntheta);
    this->spline_BR   = gsl_spline2d_alloc(gsl_interp2d_bicubic, this->npsi, this->ntheta);
    this->spline_BZ   = gsl_spline2d_alloc(gsl_interp2d_bicubic, this->npsi, this->ntheta);
    this->spline_Bphi = gsl_spline2d_alloc(gsl_interp2d_bicubic, this->npsi, this->ntheta);

    gsl_spline2d_init(this->spline_R,    this->input_r, this->theta, this->R, this->npsi, this->ntheta);
    gsl_spline2d_init(this->spline_Z,    this->input_r, this->theta, this->Z, this->npsi, this->ntheta);
    gsl_spline2d_init(this->spline_BR,   this->input_r, this->theta, this->dataBR, this->npsi, this->ntheta);
    gsl_spline2d_init(this->spline_BZ,   this->input_r, this->theta, this->dataBZ, this->npsi, this->ntheta);
    gsl_spline2d_init(this->spline_Bphi, this->input_r, this->theta, this->dataBphi, this->npsi, this->ntheta);

    this->isBuilt = true;
    return true;
}

/**
 * Calculate the Jacobian J at the specified radius and poloidal angle.
 *
 * r:     Minor radius.
 * theta: Poloidal angle.
 *
 * Optional return parameters:
 * _R:    Major radius coordinate in the given point.
 * _dRdt: Poloidal angle derivative of R.
 * _dZdt: Poloidal angle derivative of vertical coordinate Z.
 */
real_t NumericBRadialGridGenerator::JacobianAtTheta(
    const real_t r, const real_t theta,
    real_t *_R, real_t *_dRdt, real_t *_dZdt
) {
    real_t dRdr = gsl_spline2d_eval_deriv_x(this->spline_R, r, theta, this->acc_r, this->acc_theta);
    real_t dRdt = gsl_spline2d_eval_deriv_y(this->spline_R, r, theta, this->acc_r, this->acc_theta);

    real_t dZdr = gsl_spline2d_eval_deriv_x(this->spline_Z, r, theta, this->acc_r, this->acc_theta);
    real_t dZdt = gsl_spline2d_eval_deriv_y(this->spline_Z, r, theta, this->acc_r, this->acc_theta);

    real_t R    = this->Rp + gsl_spline2d_eval(this->spline_R, r, theta, this->acc_r, this->acc_theta);

    if (_R != nullptr) *_R = R;
    if (_dRdt != nullptr) *_dRdt = dRdt;
    if (_dZdt != nullptr) *_dZdt = dZdt;

    return R*abs(dRdr*dZdt - dRdt*dZdr);
}

/**
 * Calculate R/R0 at the given poloidal angle.
 *
 * r:     Minor radius.
 * theta: Poloidal angle.
 */
real_t NumericBRadialGridGenerator::ROverR0AtTheta(
    const real_t r, const real_t theta
) {
    return 1 + gsl_spline2d_eval(this->spline_R, r, theta, this->acc_r, this->acc_theta) / this->Rp;
}

/**
 * Calculate |\nabla r|^2 at the given poloidal angle.
 *
 * r:     Minor radius.
 * theta: Poloidal angle.
 */
real_t NumericBRadialGridGenerator::NablaR2AtTheta(
    const real_t r, const real_t theta
) {
    real_t R, dRdt, dZdt;
    real_t J = JacobianAtTheta(r, theta, &R, &dRdt, &dZdt);

    return R*R/(J*J) * (dRdt*dRdt + dZdt*dZdt);
}

