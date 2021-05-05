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

    this->Init(mf, frmt, ntheta_interp);
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

    this->rf_provided = new real_t[nr];
    for (len_t i = 0; i < nr+1; i++)
        this->rf_provided[i] = r_f[i];
    
    this->Init(mf, frmt, ntheta_interp);
}

/**
 * Common initializer.
 */
void NumericBRadialGridGenerator::Init(
    const std::string& mf, enum file_format frmt,
    const len_t ntheta_interp
) {
    this->isUpDownSymmetric = false;
	this->ntheta_interp = ntheta_interp;

	this->BtorGOverR0 = new real_t[GetNr()];
	this->BtorGOverR0_f = new real_t[GetNr()+1];
	this->psiPrimeRef = new real_t[GetNr()];
	this->psiPrimeRef_f = new real_t[GetNr()+1];

    this->acc_r = gsl_interp_accel_alloc();
    this->acc_theta = gsl_interp_accel_alloc();

    LoadMagneticFieldData(mf, frmt);
}

/**
 * Destructor.
 */
NumericBRadialGridGenerator::~NumericBRadialGridGenerator() {
    if (this->rf_provided != nullptr)
        delete [] this->rf_provided;

    if (this->spline_R != nullptr) {
        gsl_spline_free(this->spline_psi);
        gsl_spline2d_free(this->spline_R);
        gsl_spline2d_free(this->spline_Z);
        gsl_spline2d_free(this->spline_BR);
        gsl_spline2d_free(this->spline_BZ);
        gsl_spline2d_free(this->spline_Bphi);
    }

    gsl_interp_accel_free(this->acc_theta);
    gsl_interp_accel_free(this->acc_r);

	delete [] this->BtorGOverR0;
	delete [] this->BtorGOverR0_f;
}


/**
 * Make sure that the given poloidal angle is on the interval
 * [0, 2pi], and if not, shift so that it is.
 */
real_t NumericBRadialGridGenerator::_thetaBounded(const real_t theta) const {
    if (theta >= 0 && theta <= 2*M_PI)
        return theta;

    // if theta > 0
    //   floor(...) = N
    // if theta < 0
    //   floor(...) = -(N+1)
    return theta - 2*M_PI*std::floor(theta/(2*M_PI));
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
		this->R = this->addR0DataPoint(this->R, this->R, this->npsi, this->ntheta, 0);

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
        this->input_r[i] = this->R[i];

    // Turn 'R' and 'Z' into global coordinates (they're given
    // relative to flux surface otherwise)
    for (len_t i = 0; i < this->npsi*this->ntheta; i++) {
        this->R[i] += this->Rp;
        this->Z[i] += this->Zp;
    }

    // Calculate derived data
    this->dataB = new real_t[ntheta*npsi];
    for (len_t j = 0; j < this->ntheta; j++)
        for (len_t i = 0; i < this->npsi; i++) {
            len_t k = j*npsi+i;
            real_t Br = dataBR[k], Bz = dataBZ[k], Bp = dataBphi[k];

            // Magnetic field strength
            this->dataB[k] = sqrt(Br*Br + Bz*Bz + Bp*Bp);
        }

    delete d;
}

/**
 * Extend the given array with one element at r=0.
 */
real_t *NumericBRadialGridGenerator::addR0DataPoint(
	const real_t *x, const real_t *r, const len_t nr, const len_t ntheta,
    real_t c
) {
	real_t *arr = new real_t[(nr+1)*ntheta];
    if (std::isnan(c))
        c = x[0] - (r[0]-0.0)/(r[1]-r[0]) * (x[1] - x[0]);

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
            throw FVMException(
                "NumericBRadialGrid: Maximum r available in numeric magnetic "
                "field data is rMax = %.3f, but r = %.3f is required for radial "
                "grid. (Delta: %e)", this->input_r[this->npsi-1], rMax,
                std::abs(this->input_r[this->npsi-1]/rMax-1)
            );
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
    this->spline_psi  = gsl_spline_alloc(gsl_interp_steffen, this->npsi);

    const gsl_interp2d_type *splineType = gsl_interp2d_bicubic; //or ..._bilinear
    this->spline_R    = gsl_spline2d_alloc(splineType, this->npsi, this->ntheta);
    this->spline_Z    = gsl_spline2d_alloc(splineType, this->npsi, this->ntheta);
    this->spline_BR   = gsl_spline2d_alloc(splineType, this->npsi, this->ntheta);
    this->spline_BZ   = gsl_spline2d_alloc(splineType, this->npsi, this->ntheta);
    this->spline_Bphi = gsl_spline2d_alloc(splineType, this->npsi, this->ntheta);
    this->spline_B    = gsl_spline2d_alloc(splineType, this->npsi, this->ntheta);

    gsl_spline_init(this->spline_psi,    this->input_r, this->psi,   this->npsi);
    gsl_spline2d_init(this->spline_R,    this->input_r, this->theta, this->R, this->npsi, this->ntheta);
    gsl_spline2d_init(this->spline_Z,    this->input_r, this->theta, this->Z, this->npsi, this->ntheta);
    gsl_spline2d_init(this->spline_BR,   this->input_r, this->theta, this->dataBR, this->npsi, this->ntheta);
    gsl_spline2d_init(this->spline_BZ,   this->input_r, this->theta, this->dataBZ, this->npsi, this->ntheta);
    gsl_spline2d_init(this->spline_Bphi, this->input_r, this->theta, this->dataBphi, this->npsi, this->ntheta);
    gsl_spline2d_init(this->spline_B,    this->input_r, this->theta, this->dataB, this->npsi, this->ntheta);

	// Reference quantities
	for (len_t i = 0; i < GetNr(); i++) {
		real_t
			R  = gsl_spline2d_eval(
					this->spline_R, this->r[i], 0,
					this->acc_r, this->acc_theta
				),
			Bphi = gsl_spline2d_eval(
					this->spline_Bphi, this->r[i], 0,
					this->acc_r, this->acc_theta
				),
            psip = gsl_spline_eval_deriv(
                    this->spline_psi, this->r[i], this->acc_r
                );

		this->BtorGOverR0[i] = (R/this->Rp) * Bphi;
        this->psiPrimeRef[i] = psip;
		//this->psiPrimeRef[i] = (R/this->Rp) * hypot(Br, Bz);
	}
	for (len_t i = 0; i < GetNr()+1; i++) {
		real_t
			R  = gsl_spline2d_eval(
					this->spline_R, this->r_f[i], 0,
					this->acc_r, this->acc_theta
				),
			Bphi = gsl_spline2d_eval(
					this->spline_Bphi, this->r_f[i], 0,
					this->acc_r, this->acc_theta
				),
            psip = gsl_spline_eval_deriv(
                    this->spline_psi, this->r_f[i], this->acc_r
                );

		this->BtorGOverR0_f[i] = (R/this->Rp) * Bphi;
        this->psiPrimeRef_f[i] = psip;
		//this->psiPrimeRef_f[i] = (R/this->Rp) * hypot(Br, Bz);
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
 * Calculate the Jacobian J, normalized to R0, at the specified radius and
 * poloidal angle.
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
    real_t t = this->_thetaBounded(theta);

    real_t dRdr = gsl_spline2d_eval_deriv_x(this->spline_R, r, t, this->acc_r, this->acc_theta);
    real_t dRdt = gsl_spline2d_eval_deriv_y(this->spline_R, r, t, this->acc_r, this->acc_theta);

    real_t dZdr = gsl_spline2d_eval_deriv_x(this->spline_Z, r, t, this->acc_r, this->acc_theta);
    real_t dZdt = gsl_spline2d_eval_deriv_y(this->spline_Z, r, t, this->acc_r, this->acc_theta);

    real_t R    = gsl_spline2d_eval(this->spline_R, r, t, this->acc_r, this->acc_theta);

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
    real_t t = this->_thetaBounded(theta);
    return gsl_spline2d_eval(this->spline_R, r, t, this->acc_r, this->acc_theta) / this->Rp;
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
    // Avoid division by zero
    // (this value is technically incorrect, but since NablaR2AtTheta
    // at r=0 is never used anywhere, we do this to avoid the hassle
    // of extrapolating...)
	if (r == 0)
		return 1.0;

    real_t R, dRdt, dZdt;
    real_t J_R0 = JacobianAtTheta(r, theta, &R, &dRdt, &dZdt);

	return R*R/(Rp*Rp*J_R0*J_R0) * (dRdt*dRdt + dZdt*dZdt);
}

/**
 * Evaluate all the geometric quantities in one go.
 */
void NumericBRadialGridGenerator::EvaluateGeometricQuantities(
    const real_t r, const real_t theta,
    real_t &B, real_t &Jacobian, real_t &ROverR0,
    real_t &NablaR2
) {
    real_t t = this->_thetaBounded(theta);
    real_t R, dRdt, dZdt;

    Jacobian = JacobianAtTheta(r, theta, &R, &dRdt, &dZdt);
    ROverR0  = R/this->Rp;
    NablaR2  = r==0 ? 0 : (R*R/(Jacobian*Jacobian) * (dRdt*dRdt + dZdt*dZdt));

    B = gsl_spline2d_eval(this->spline_B, r, t, this->acc_r, this->acc_theta);
}

/**
 * Evaluate magnetic field strength B at given poloidal angle
 * and radius.
 */
real_t NumericBRadialGridGenerator::EvalB(const real_t r, const real_t theta) {
    real_t t    = this->_thetaBounded(theta);
    return gsl_spline2d_eval(this->spline_B, r, t, this->acc_r, this->acc_theta);
}

real_t NumericBRadialGridGenerator::BAtTheta(const len_t ir, const real_t theta) {
	return EvalB(this->r[ir], theta);
}
real_t NumericBRadialGridGenerator::BAtTheta_f(const len_t ir, const real_t theta) {
	return EvalB(this->r_f[ir], theta);
}

/**
 * Calculate minor radius coordinate 'r' corresponding to the given
 * Cartesian coordinates (x,y,z).
 *
 * (The Cartesian coordinate system is oriented such that x and y span
 * the poloidal plane. The origin of x and y is the magnetic axis.)
 */
void NumericBRadialGridGenerator::GetRThetaFromCartesian(real_t *r, real_t *theta,
    real_t x, real_t y, real_t z, real_t lengthScale
) {
    // Major radius coordinate
    real_t  R = hypot(x-R0, z);

    // Position vector
    real_t rhox = x-R0 - R0*(x-R0)/R;
    real_t rhoy = y;
    real_t rhoz = z-R0 - R0*(z-R0)/R;

    // Minor radius at poloidal angle
    real_t rho = sqrt(rhox*rhox + rhoy*rhoy + rhoz*rhoz);

    // Poloidal angle
    if (R >= R0)
        *theta = std::atan2(rhoy, +hypot(rhox, rhoz));
    else
        *theta = std::atan2(rhoy, -hypot(rhox, rhoz));

    // Bisection to find radial coordinate corresponding
    // to 'r' at 'theta'...
    int nr=GetNr();
    real_t ra = 0, rb=this->r_f[nr-1];
    do {
        *r = (ra-rb)/2;
        real_t
            xx = gsl_spline2d_eval(
                this->spline_R, *r, *theta,
                this->acc_r, this->acc_theta
            ),
            yy = gsl_spline2d_eval(
                this->spline_Z, *r, *theta,
                this->acc_r, this->acc_theta
            );

        if (hypot(xx, yy) < rho)
            ra = *r;
        else
            rb = *r;
    } while(std::abs(rb-ra) > lengthScale*1e-3);
}

/**
 * ???
 */
void NumericBRadialGridGenerator::GetGradRCartesian(real_t*, real_t, real_t) {
}

/**
 * ???
 */
real_t NumericBRadialGridGenerator::FindClosestApproach(
    real_t, real_t, real_t,
    real_t, real_t, real_t
) {
    return 0;
}
