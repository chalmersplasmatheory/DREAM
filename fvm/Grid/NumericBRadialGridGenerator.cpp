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
 * nr:            Number of radial grid points in uniform radial *distribution* grid.
 * r0:            Value of innermost point on radial *flux* grid.
 * ra:            Value of outermost point on radial *flux* grid.
 * mf:            Name of file containing the magnetic field data to load.
 * frmt:          Format in which the magnetic field data is stored.
 * ntheta_interp: Poloidal angle resolution in quadrature for flux surface and bounce averages.
 */
NumericBRadialGridGenerator::NumericBRadialGridGenerator(
    const len_t nr, const real_t r0, const real_t ra,
    const std::string& mf, enum file_format frmt,
	const len_t ntheta_interp
) : RadialGridGenerator(nr), rMin(r0), rMax(ra) {

	this->ntheta_interp = ntheta_interp;

	this->BtorGOverR0 = new real_t[nr];
	this->BtorGOverR0_f = new real_t[nr+1];
	this->psiPrimeRef = new real_t[nr];
	this->psiPrimeRef_f = new real_t[nr+1];

    this->isUpDownSymmetric = false;
    LoadMagneticFieldData(mf, frmt);
}

/**
 * Constructor for a custom radial grid.
 *
 * r_f:           Radial flux grid desired.
 * nr:            Number of points on distribution grids corresponding to 
 *                'r_f' (i.e. the number of points in 'r_f', minus one).
 * ntheta_interp: Poloidal angle resolution in quadrature for flux surface and bounce averages.
 */
NumericBRadialGridGenerator::NumericBRadialGridGenerator(
    const real_t *r_f, const len_t nr,
    const std::string& mf, enum file_format frmt,
	const len_t ntheta_interp
) : RadialGridGenerator(nr) {

    this->isUpDownSymmetric = false;
	this->ntheta_interp = ntheta_interp;

    this->rf_provided = new real_t[nr];
    for (len_t i = 0; i < nr+1; i++)
        this->rf_provided[i] = r_f[i];

	this->BtorGOverR0 = new real_t[nr];
	this->BtorGOverR0_f = new real_t[nr+1];
	this->psiPrimeRef = new real_t[nr];
	this->psiPrimeRef_f = new real_t[nr+1];

    LoadMagneticFieldData(mf, frmt);
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

	delete [] this->BtorGOverR0;
	delete [] this->BtorGOverR0_f;
}


/**
 * Load magnetic field data from the named file.
 *
 * filename: Name of file to load data from.
 */
void NumericBRadialGridGenerator::LoadMagneticFieldData(
    const std::string& filename, enum file_format frmt
) {
    SFile *sf = SFile::Create(filename, SFILE_MODE_READ);
    this->LoadMagneticFieldData(sf, frmt);

    sf->Close();
    delete sf;
}

/**
 * Load magnetic field data from the given file.
 *
 * sf: SFile object representing the file.
 */
void NumericBRadialGridGenerator::LoadMagneticFieldData(
    SFile *sf, enum file_format frmt
) {
    struct NumericBData *d;
    switch (frmt) {
        case FILE_FORMAT_LUKE:
            d = LoadNumericBFromLUKE(sf);
            break;

        default:
            throw FVMException(
                "NumericBRadialGrid: Unrecognized file format specified for magnetic field data: %d.",
                frmt
            );
    }

    auto convert_data = [](const double *a, const len_t nx, const len_t ny) {
        real_t *d = new real_t[nx*ny];
        for (len_t i = 0; i < nx*ny; i++)
            d[i] = (real_t)a[i];

        return d;
    };

	this->npsi = d->npsi;
	this->ntheta = d->ntheta;

	// Set plasma parameters
	this->Rp = d->Rp;
	this->Zp = d->Zp;

    // Set magnetic field data
    this->psi      = convert_data(d->psi, npsi, 1);
    this->theta    = convert_data(d->theta, ntheta, 1);
    
    this->R        = convert_data(d->R, ntheta, npsi);
    this->Z        = convert_data(d->Z, ntheta, npsi);

    this->dataBR   = convert_data(d->Br, ntheta, npsi);
    this->dataBZ   = convert_data(d->Bz, ntheta, npsi);
    this->dataBphi = convert_data(d->Bphi, ntheta, npsi);

	// Add point at r=0 if it does not already exist...
	if (this->R[0] != 0) {
		this->psi = this->addR0DataPoint(this->psi, this->R, this->npsi, 1);
		this->dataBR = this->addR0DataPoint(this->dataBR, this->R, this->npsi, this->ntheta);
		this->dataBZ = this->addR0DataPoint(this->dataBZ, this->R, this->npsi, this->ntheta);
		this->dataBphi = this->addR0DataPoint(this->dataBphi, this->R, this->npsi, this->ntheta);

		this->Z = this->addR0DataPoint(this->Z, this->R, this->npsi, this->ntheta);
		// The radial grid should be modified last...
		this->R = this->addR0DataPoint(this->R, this->R, this->npsi, this->ntheta);

		// r grid has now been extended...
		this->npsi++;
	}

	// Ensure that theta=0 exists...
	if (this->theta[0] != 0)
		throw FVMException("NumericBRadialGrid: All numerical data must be given in theta = 0.");

	// Add point at theta=2*pi if it does not already exist...
	if (this->theta[this->ntheta-1] != 2*M_PI) {
		this->theta = this->addThetaDataPoint(this->theta, 1, this->ntheta);
		this->theta[this->ntheta] = 2*M_PI;
		this->R = this->addThetaDataPoint(this->R, this->npsi, this->ntheta);
		this->Z = this->addThetaDataPoint(this->Z, this->npsi, this->ntheta);
		this->dataBR = this->addThetaDataPoint(this->dataBR, this->npsi, this->ntheta);
		this->dataBZ = this->addThetaDataPoint(this->dataBZ, this->npsi, this->ntheta);
		this->dataBphi = this->addThetaDataPoint(this->dataBphi, this->npsi, this->ntheta);

		// r grid has now been extended...
		this->ntheta++;
	}

    // Evaluate minor radius in outer midplane
    this->input_r = new real_t[this->npsi];
    for (len_t i = 0; i < this->npsi; i++)
        this->input_r[i] = this->R[i] - this->Rp;

    delete d;
}

/**
 * Extend the given array with one element at r=0.
 */
real_t *NumericBRadialGridGenerator::addR0DataPoint(
	const real_t *x, const real_t *r, const len_t nr, const len_t ntheta
) {
	real_t *arr = new real_t[(nr+1)*ntheta];
	real_t c = x[0] - (r[0]-this->Rp)/(r[1]-r[0]) * (x[1] - x[0]);

	for (len_t j = 0; j < ntheta; j++) {
		// Value at r=0
		arr[j*(nr+1)] = c;

		for (len_t i = 0; i < nr; i++)
			arr[j*(nr+1) + i+1] = x[j*nr+i];
	}

	delete [] x;

	return arr;
}

/**
 * Extend the given array with one element at theta=2*pi.
 */
real_t *NumericBRadialGridGenerator::addThetaDataPoint(
	const real_t *x, const len_t nr, const len_t ntheta
) {
	real_t *arr = new real_t[nr*(ntheta+1)];

	for (len_t i = 0; i < ntheta*nr; i++)
		arr[i] = x[i];
	
	// Add x(theta=2*pi) = x(theta=0)
	for (len_t i = 0; i < nr; i++)
		arr[ntheta*nr + i] = x[i];
	
	delete [] x;

	return arr;
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
            throw FVMException("NumericBRadialGrid: Maximum r available in numeric magnetic field data is rMax = %.3f, but r = %.3f is required for radial grid.", this->input_r[this->npsi-1], rMax);
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
    this->spline_BR   = gsl_spline2d_alloc(gsl_interp2d_bilinear, this->npsi, this->ntheta);
    this->spline_BZ   = gsl_spline2d_alloc(gsl_interp2d_bilinear, this->npsi, this->ntheta);
    this->spline_Bphi = gsl_spline2d_alloc(gsl_interp2d_bilinear, this->npsi, this->ntheta);

    gsl_spline2d_init(this->spline_R,    this->input_r, this->theta, this->R, this->npsi, this->ntheta);
    gsl_spline2d_init(this->spline_Z,    this->input_r, this->theta, this->Z, this->npsi, this->ntheta);
    gsl_spline2d_init(this->spline_BR,   this->input_r, this->theta, this->dataBR, this->npsi, this->ntheta);
    gsl_spline2d_init(this->spline_BZ,   this->input_r, this->theta, this->dataBZ, this->npsi, this->ntheta);
    gsl_spline2d_init(this->spline_Bphi, this->input_r, this->theta, this->dataBphi, this->npsi, this->ntheta);

	// Reference quantities
	for (len_t i = 0; i < GetNr(); i++) {
		real_t
			R  = gsl_spline2d_eval(
					this->spline_R, this->r[i], 0,
					this->acc_r, this->acc_theta
				),
			Br = gsl_spline2d_eval(
					this->spline_BR, this->r[i], 0,
					this->acc_r, this->acc_theta
				),
			Bz = gsl_spline2d_eval(
					this->spline_BZ, this->r[i], 0,
					this->acc_r, this->acc_theta
				),
			Bphi = gsl_spline2d_eval(
					this->spline_Bphi, this->r[i], 0,
					this->acc_r, this->acc_theta
				);

		this->BtorGOverR0[i] = (R/this->Rp) * Bphi;
		this->psiPrimeRef[i] = (R/this->Rp) * hypot(Br, Bz);
	}
	for (len_t i = 0; i < GetNr()+1; i++) {
		real_t
			R  = gsl_spline2d_eval(
					this->spline_R, this->r_f[i], 0,
					this->acc_r, this->acc_theta
				),
			Br = gsl_spline2d_eval(
					this->spline_BR, this->r_f[i], 0,
					this->acc_r, this->acc_theta
				),
			Bz = gsl_spline2d_eval(
					this->spline_BZ, this->r_f[i], 0,
					this->acc_r, this->acc_theta
				),
			Bphi = gsl_spline2d_eval(
					this->spline_Bphi, this->r_f[i], 0,
					this->acc_r, this->acc_theta
				);

		this->BtorGOverR0_f[i] = (R/this->Rp) * Bphi;
		this->psiPrimeRef_f[i] = (R/this->Rp) * hypot(Br, Bz);
	}
	
	rGrid->SetReferenceMagneticFieldData(
		this->BtorGOverR0, this->BtorGOverR0_f,
		this->psiPrimeRef, this->psiPrimeRef_f,
		this->Rp
	);

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

    return R/this->Rp*fabs(dRdr*dZdt - dRdt*dZdr);
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
	if (r == 0)
		return 0.0;

    real_t R, dRdt, dZdt;
    real_t J = JacobianAtTheta(r, theta, &R, &dRdt, &dZdt);

	return R*R/(J*J) * (dRdt*dRdt + dZdt*dZdt);
}

/**
 * Evaluate all the geometric quantities in one go.
 */
void NumericBRadialGridGenerator::EvaluateGeometricQuantities(
    const real_t r, const real_t theta,
    real_t &B, real_t &Jacobian, real_t &ROverR0,
    real_t &NablaR2
) {
    real_t R, dRdt, dZdt;

    Jacobian = JacobianAtTheta(r, theta, &R, &dRdt, &dZdt);
    ROverR0  = R/this->Rp;
    NablaR2  = r==0 ? 0 : (R*R/(Jacobian*Jacobian) * (dRdt*dRdt + dZdt*dZdt));

    real_t Br, Bz, Bphi;
    Br   = gsl_spline2d_eval(this->spline_BR, r, theta, this->acc_r, this->acc_theta);
    Bz   = gsl_spline2d_eval(this->spline_BZ, r, theta, this->acc_r, this->acc_theta);
    Bphi = gsl_spline2d_eval(this->spline_Bphi, r, theta, this->acc_r, this->acc_theta);

    B    = sqrt(Br*Br + Bz*Bz + Bphi*Bphi);
}

/**
 * Evaluate magnetic field strength B at given poloidal angle
 * and radius.
 */
real_t NumericBRadialGridGenerator::EvalB(const real_t r, const real_t theta) {
    real_t Br   = gsl_spline2d_eval(this->spline_BR, r, theta, this->acc_r, this->acc_theta);
    real_t Bz   = gsl_spline2d_eval(this->spline_BZ, r, theta, this->acc_r, this->acc_theta);
    real_t Bphi = gsl_spline2d_eval(this->spline_Bphi, r, theta, this->acc_r, this->acc_theta);

    return sqrt(Br*Br + Bz*Bz + Bphi*Bphi);
}

real_t NumericBRadialGridGenerator::BAtTheta(const len_t ir, const real_t theta) {
	return EvalB(this->r[ir], theta);
}
real_t NumericBRadialGridGenerator::BAtTheta_f(const len_t ir, const real_t theta) {
	return EvalB(this->r_f[ir], theta);
}

/**
 * Locate a magnetic field extremum point.
 *
 * ir:    Radial grid index to locate extremum at.
 * sgn:   +1: locate minimum, -1: locate maximum.
 * fgrid: Type of grid to look for extremum on (distribution or flux grid).
 *
 * NOTE: This routine assumes that the magnetic field only has global
 * extrema on all flux surfaces.
 */
real_t NumericBRadialGridGenerator::FindMagneticFieldExtremum(
	len_t ir, int_t sgn, enum fluxGridType fgrid
) {
	const real_t EPSABS = 1e-6;
	real_t lower, upper, guess;

	// Determine if we are going to search on [0,pi] or [pi,2*pi]
	real_t B, Beps;
	if (sgn > 0)
		guess = 0;
	else
		guess = M_PI;

	// Evaluate magnetic field
	if (fgrid == FLUXGRIDTYPE_DISTRIBUTION) {
		B = BAtTheta(ir, guess+sgn*EPSABS);
		Beps = BAtTheta(ir, 2*M_PI-(guess+sgn*EPSABS));
	} else {
		B = BAtTheta_f(ir, guess+sgn*EPSABS);
		Beps = BAtTheta_f(ir, 2*M_PI-(guess+sgn*EPSABS));
	}

	// Is extremum in upper half plane?
	if (sgn*B < sgn*Beps)
		lower = 0, upper = M_PI;
	else
		lower = M_PI, upper = 2*M_PI;

	return FindMagneticFieldExtremum_inner(ir, sgn, fgrid, lower, upper);
}

