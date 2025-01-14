/**
 * Implementation of the numeric magnetic field radial grid generator. This
 * grid generator loads a numeric magnetic field from the specified file and
 * builds a correspondingly shaped radial grid.
 */
#include <algorithm>
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

    this->rf_provided = new real_t[nr+1];
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
	
	if (this->input_r != nullptr)
		delete [] this->input_r;
	
	if (this->dataB != nullptr)
		delete [] this->dataB;
	
	if (this->theta != nullptr)
		delete [] this->theta;
	if (this->R != nullptr)
		delete [] this->R;
	if (this->Z != nullptr)
		delete [] this->Z;
	if (this->psi != nullptr)
		delete [] this->psi;
	if (this->dataBR != nullptr)
		delete [] this->dataBR;
	if (this->dataBZ != nullptr)
		delete [] this->dataBZ;
	if (this->dataBphi != nullptr)
		delete [] this->dataBphi;

    if (this->spline_R != nullptr) {
        gsl_spline_free(this->spline_psi);
        gsl_spline2d_free(this->spline_R);
        gsl_spline2d_free(this->spline_Z);
        gsl_spline2d_free(this->spline_BR);
        gsl_spline2d_free(this->spline_BZ);
        gsl_spline2d_free(this->spline_Bphi);
		gsl_spline2d_free(this->spline_B);
    }

    gsl_interp_accel_free(this->acc_theta);
    gsl_interp_accel_free(this->acc_r);
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

		this->Z = this->addR0DataPoint(this->Z, this->R, this->npsi, this->ntheta, 0); // <----
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
	} else {
		// Ensure that values at theta=0 and theta=2*pi are identical
		for (len_t ir = 0; ir < this->npsi; ir++) {
			len_t idx0 = 0*this->npsi + ir;
			len_t idx1 = (this->ntheta-1)*this->npsi + ir;

			if (this->R[idx0] != this->R[idx1])
				this->R[idx1] = this->R[idx0];
			if (this->Z[idx0] != this->Z[idx1])
				this->Z[idx1] = this->Z[idx0];
			if (this->dataBR[idx0] != this->dataBR[idx1])
				this->dataBR[idx1] = this->dataBR[idx0];
			if (this->dataBZ[idx0] != this->dataBZ[idx1])
				this->dataBZ[idx1] = this->dataBZ[idx0];
			if (this->dataBphi[idx0] != this->dataBphi[idx1])
				this->dataBphi[idx1] = this->dataBphi[idx0];
		}
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

    // Verify that the magnetic field contains exactly one maximum and
    // one minimum

    // Tolelance for an extremum to be announced...
    const real_t TOLERANCE = sqrt(std::numeric_limits<real_t>::epsilon());
    for (len_t i = 0; i < this->npsi; i++) {
        bool minFound = false, maxFound = false;

		len_t maxk=0, mink=0;
        for (len_t j = 1; j < this->ntheta; j++) {
            len_t k  = j*npsi+i;
            len_t km = k-npsi, kp = k+npsi;

            // Wrap-around at the upper endpoint
            if (j == this->ntheta-1)
                kp = 0*npsi + i;

            real_t dB =
                std::max(
                    std::abs(this->dataB[k]-this->dataB[km]),
                    std::abs(this->dataB[kp]-this->dataB[k])
                );

            // Is this a (local) maximum or minimum?
            if (this->dataB[km] < this->dataB[k] &&
                this->dataB[kp] < this->dataB[k] &&
                dB > TOLERANCE) {

                // Has a maximum already been found?
                if (maxFound)
                    throw FVMException(
                        "The numeric magnetic field has more than one maximum "
                        "along at least one magnetic field line."
						"ipsi = " LEN_T_PRINTF_FMT
						", itheta(1) = " LEN_T_PRINTF_FMT
						", itheta(2) = " LEN_T_PRINTF_FMT,
						i, (k-i)/npsi, (maxk-i)/npsi
                    );
                else {
                    maxFound = true;
					maxk = k;
				}
            } else if (this->dataB[km] > this->dataB[k] &&
                       this->dataB[kp] > this->dataB[k] &&
                       dB > TOLERANCE) {
                
                // Has a minimum already been found?
                if (minFound)
                    throw FVMException(
                        "The numeric magnetic field has more than one minimum "
                        "along at least one magnetic field line. "
						"ipsi = " LEN_T_PRINTF_FMT
						", itheta(1) = " LEN_T_PRINTF_FMT
						", itheta(2) = " LEN_T_PRINTF_FMT,
						i, (k-i)/npsi, (mink-i)/npsi
                    );
                else {
                    minFound = true;
					mink = k;
				}
            }
        }
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
        if (rf_provided[0] < 0)
            throw FVMException("NumericBRadialGrid: First point on custom radial grid is less than zero.");
        else if (rf_provided[GetNr()] > this->input_r[this->npsi-1])
            throw FVMException(
                "NumericBRadialGrid: Last point on custom radial grid may not be greater than "
                "the maximum r available in the numeric magnetic field data, rMax = %.3f.",
                this->input_r[this->npsi-1]
            );
        else if (rf_provided[0] >= rf_provided[GetNr()])
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
void NumericBRadialGridGenerator::GetRThetaPhiFromCartesian(real_t *r, real_t *theta, real_t *phi,
    real_t x, real_t y, real_t z, real_t lengthScale, real_t startingGuessR
) {

    // If x+Rp<0, the cartesian coordinates are clearly outside the radial grid
    // and we therefore immediately return an arbitrary value outside the radial grid
    if (x+Rp<0){
        *theta = 0;
        *phi = 0;
        *r=this->r_f[GetNr()]+1e-2;
        return;
    }
    // Major radius coordinate
    real_t  R = hypot(x+Rp, z);

    // Position vector
    real_t rhox = x+Rp - Rp*(x+Rp)/R;
    real_t rhoy = y;
    real_t rhoz = z - Rp*z/R;

    // Minor radius at poloidal angle
    real_t rho = sqrt(rhox*rhox + rhoy*rhoy + rhoz*rhoz);
    
    real_t r_tmp;
    real_t theta_tmp;

    // Poloidal angle
    if (R >= Rp)
        theta_tmp = std::atan2(rhoy, +hypot(rhox, rhoz));
    else
        theta_tmp = std::atan2(rhoy, -hypot(rhox, rhoz));
        
    if (theta_tmp < 0)
        theta_tmp+=2*M_PI;
        
    *theta=theta_tmp;
        

	// Bisection to find radial coordinate corresponding
	// to 'r' at 'theta'...
	// We make a guess for a valid search intervall of startingGuessR+/-lengthScale, 
	// and check if it has to be expanded before actually starting with the bisection
	real_t ra = std::max(0.0, startingGuessR-lengthScale);
	real_t rb = ra + 2*lengthScale;
	if(rb>r_f[GetNr()]){
	    ra = std::max(0.0, ra -( rb - r_f[GetNr()])); 
	    rb = r_f[GetNr()]; 
	}
	real_t rhoa, rhob;
	do {
		real_t
		    xxa = gsl_spline2d_eval(
		        this->spline_R, ra, *theta,
		        this->acc_r, this->acc_theta
		    ) - Rp,
		    yya = gsl_spline2d_eval(
		        this->spline_Z, ra, *theta,
		        this->acc_r, this->acc_theta
		    ) - this->Zp;
		   
		real_t
		    xxb = gsl_spline2d_eval(
		        this->spline_R, rb, *theta,
		        this->acc_r, this->acc_theta
		    )-Rp,
		    yyb = gsl_spline2d_eval(
		        this->spline_Z, rb, *theta,
		        this->acc_r, this->acc_theta
		    )- this->Zp;
		rhoa=hypot(xxa,yya);
		rhob=hypot(xxb,yyb);
		if(rhoa>rho && rhob>rho){
			ra=std::max(0.0, ra-2*lengthScale);
			rb=ra + 2*lengthScale;
		}
	    else if(rhoa<rho && rhob<rho){
	        ra+=2*lengthScale;
	        rb+=2*lengthScale;
	    }
	    if(rb>this->r_f[GetNr()])
	        break;
	  } while ((rhoa>rho && rhob>rho) || (rhoa<rho && rhob<rho));
	  
	// Make the bisection
	if(rb<this->r_f[GetNr()]){
	    do {
	        r_tmp = (ra+rb)/2;
	        real_t
	            xx = gsl_spline2d_eval(
	                this->spline_R, r_tmp, *theta,
	                this->acc_r, this->acc_theta
	            )-Rp,
	            yy = gsl_spline2d_eval(
	                this->spline_Z, r_tmp, *theta,
	                this->acc_r, this->acc_theta
	            ) - this->Zp;

	        if (hypot(xx, yy) < rho)
	            ra = r_tmp;
	        else
	            rb = r_tmp;
	    } while(std::abs(rb-ra) > lengthScale*CartesianCoordinateTol);
	    
        *r=r_tmp;
    } else{
        *r=this->r_f[GetNr()]+1e-2; // Arbitrary value outside the radial grid
    }

	// Newton solver (disabled for now)
	/*r_tmp=startingGuessR;
	do {
	    real_t
	        xx = gsl_spline2d_eval(
	            this->spline_R, r_tmp, *theta,
	            this->acc_r, this->acc_theta
	        ),
	        yy = gsl_spline2d_eval(
	            this->spline_Z, r_tmp, *theta,
	            this->acc_r, this->acc_theta
	        );
	        
	    // note that r_tmp is here the x-variable given to gsl_spline2d_eval_deriv_x, 
	    // so that it should really be the x-derivative which is evaluated for both xx and yy
	    real_t
	        dxxdr = gsl_spline2d_eval_deriv_x(
	            this->spline_R, r_tmp, *theta,
	            this->acc_r, this->acc_theta
	        ),
	        dyydr = gsl_spline2d_eval_deriv_x(
	            this->spline_Z, r_tmp, *theta,
	            this->acc_r, this->acc_theta
	        );
	        
	    real_t rho_newton=hypot(xx,yy);
	    r_tmp=r_tmp-(rho_newton-rho)/((xx*dxxdr+yy*dyydr)/rho_newton);
	} while(std::abs(rho_newton-rho) > lengthScale * tolFactor);
	*r=r_tmp;*/
	
	*phi = atan2(z,(Rp+x)); 
	
}

/**
 * Calculates the gradient of the minor radius coordinate 'r' in cartesian coordinates
 */
void NumericBRadialGridGenerator::GetGradRCartesian(real_t* gradr, real_t r, real_t theta, real_t phi) {
	//throw FVMException("NumericBRadialGridGenerator: This module is currently incompatible with the SPI module.");
	if(r<r_f[GetNr()]){
        real_t
        dRdr = gsl_spline2d_eval_deriv_x(
            this->spline_R, r, theta,
            this->acc_r, this->acc_theta
        ),
        dzdr = gsl_spline2d_eval_deriv_x(
            this->spline_Z, r, theta,
            this->acc_r, this->acc_theta
        ),    
        dRdtheta = gsl_spline2d_eval_deriv_y(
            this->spline_R, r, theta,
            this->acc_r, this->acc_theta
        ),    
        dzdtheta = gsl_spline2d_eval_deriv_y(
            this->spline_Z, r, theta,
            this->acc_r, this->acc_theta
        );
        
        real_t common_factor = 1/(dRdr*dzdtheta - dRdtheta*dzdr);
        gradr[0] = common_factor * dzdtheta * cos(phi);
        gradr[1] = - common_factor * dRdtheta;
        gradr[2] = common_factor * dzdtheta * sin(phi);
    }else{
        gradr[0] = 0;
        gradr[1] = 0;
        gradr[2] = 0;
    }
}


/**
 * Return a list of flux surface R coordinates
 * on the simulation radial grid.
 */
const real_t *NumericBRadialGridGenerator::GetFluxSurfaceRMinusR0() {
	const len_t nr = this->GetNr();
	real_t *R = new real_t[nr * this->ntheta];

	for (len_t j = 0, i = 0; j < ntheta; j++)
		for (len_t ir = 0; ir < nr; ir++, i++)
			R[i] = ROverR0AtTheta(ir, this->theta[j]) * this->Rp - this->Rp;

	return R;
}


/**
 * Return a list of flux surface R coordinates
 * on the simulation radial grid.
 */
const real_t *NumericBRadialGridGenerator::GetFluxSurfaceRMinusR0_f() {
	const len_t nr = this->GetNr();
	real_t *R = new real_t[(nr+1) * this->ntheta];

	for (len_t j = 0, i = 0; j < ntheta; j++)
		for (len_t ir = 0; ir < nr+1; ir++, i++)
			R[i] = ROverR0AtTheta_f(ir, this->theta[j]) * this->Rp - this->Rp;

	return R;
}


/**
 * Returns a list of flux surface Z coordinates
 * on the simulation grid.
 */
const real_t *NumericBRadialGridGenerator::GetFluxSurfaceZMinusZ0() {
	const len_t nr = this->GetNr();
	real_t *Z = new real_t[nr * this->ntheta];

	for (len_t j = 0, i = 0; j < ntheta; j++) {
		for (len_t ir = 0; ir < nr; ir++, i++) {
			Z[i] = gsl_spline2d_eval(
				this->spline_Z, this->r[ir], this->theta[j],
				this->acc_r, this->acc_theta
			) - this->Zp;
		}
	}

	return Z;
}


/**
 * Returns a list of flux surface Z coordinates
 * on the simulation grid.
 */
const real_t *NumericBRadialGridGenerator::GetFluxSurfaceZMinusZ0_f() {
	const len_t nr = this->GetNr();
	real_t *Z = new real_t[(nr+1) * this->ntheta];

	for (len_t j = 0, i = 0; j < ntheta; j++) {
		for (len_t ir = 0; ir < nr+1; ir++, i++) {
			Z[i] = gsl_spline2d_eval(
				this->spline_Z, this->r_f[ir], this->theta[j],
				this->acc_r, this->acc_theta
			) - this->Zp;
		}
	}

	return Z;
}


/**
 * Returns a list of poloidal angles on which the flux
 * surfaces returned by 'GetFluxSurfaceR()' and 'GetFluxSurfaceZ()'
 * are defined.
 */
const real_t *NumericBRadialGridGenerator::GetPoloidalAngle() {
	real_t *theta = new real_t[ntheta];
	for (len_t i = 0; i < ntheta; i++)
		theta[i] = this->theta[i];
	
	return theta;
}


/*
 * Save the magnitude of the magnetic field vector to the named
 * output file (saved using the 'SFile' API).
 *
 * filename: Name of file to save data to.
 */
void NumericBRadialGridGenerator::__SaveB(const char *filename) {
    // DEBUG: Save magnetic field
    const len_t NTHETA = 1000;
    real_t **B = new real_t*[GetNr()];
    B[0] = new real_t[GetNr()*NTHETA];
    for (len_t i = 0; i < GetNr(); i++) {
        if (i > 0)
            B[i] = B[i-1] + NTHETA;

        for (len_t j = 0; j < NTHETA; j++)
            B[i][j] = this->BAtTheta(i, j*2*M_PI/NTHETA);
    }

    SFile *sf = SFile::Create(filename, SFILE_MODE_WRITE);
    sf->WriteArray("B", B, GetNr(), NTHETA);
    sf->Close();
}

