/**
 * Implementation of the numeric magnetic field radial grid generator. This
 * grid generator loads a numeric magnetic field from the specified file and
 * builds a correspondingly shaped radial grid.
 */
#include <algorithm>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline2d.h>
#include "FVM/Grid/Stellarator/NumericBRadialGridGenerator.hpp"


using namespace DREAM::FVM;


/**
 * Constructor for a uniform radial grid.
 *
 * nr:            Number of radial grid points in uniform radial *distribution* grid. 
 * r0:            Value of innermost point on radial *flux* grid.                     // TODO: Remove?
 * ra:            Value of outermost point on radial *flux* grid.                     // TODO: Remove?
 * mf:            Name of file containing the magnetic field and Jacobian data to load.
 * frmt:          Format in which the magnetic field data is stored.                  // TODO: DESC
 * ntheta_interp: Poloidal angle resolution in quadrature for flux surface and bounce averages.
 * nphi_interp:   Toroidal angle resolution in quadrature for flux surface and bounce averages.
 */
NumericBRadialGridGenerator::NumericBRadialGridGenerator(
    const len_t nr, const real_t r0, const real_t ra,
    const std::string& mf, enum file_format frmt,
	const len_t ntheta_interp, const len_t nphi_interp
) : RadialGridGenerator(nr), rMin(r0), rMax(ra) {
    const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2; // TODO: Ok?
    gsl_multi_fmin = gsl_multimin_fminimizer_alloc(T, 2);

    this->Init(mf, frmt, ntheta_interp, nphi_interp);
}

/**
 * Constructor for a custom radial grid.
 *
 * r_f:           Radial flux grid desired.                                           // TODO: Adapt?
 * nr:            Number of points on distribution grids corresponding to 
 *                'r_f' (i.e. the number of points in 'r_f', minus one).
 * mf:            Name of file containing the magnetic field and Jacobian data to load.
 * frmt:          Format in which the magnetic field data is stored.                  // TODO: DESC
 * ntheta_interp: Poloidal angle resolution in quadrature for flux surface and bounce averages.
 * nphi_interp:   Toroidal angle resolution in quadrature for flux surface and bounce averages.
 */
NumericBRadialGridGenerator::NumericBRadialGridGenerator(
    const real_t *r_f, const len_t nr,
    const std::string& mf, enum file_format frmt,
	const len_t ntheta_interp, const len_t nphi_interp
) : RadialGridGenerator(nr) {

    this->rf_provided = new real_t[nr+1];
    for (len_t i = 0; i < nr+1; i++)
        this->rf_provided[i] = r_f[i];
    
    this->Init(mf, frmt, ntheta_interp, nphi_interp);
}

/**
 * Common initializer.
 */
void NumericBRadialGridGenerator::Init(
    const std::string& mf, enum file_format frmt,
    const len_t ntheta_interp, const len_t nphi_interp
) {
    this->nfp = 0; // TODO
	this->ntheta_interp = ntheta_interp;
	this->nphi_interp   = nphi_interp;

	this->BtorGOverR0 = new real_t[GetNr()];
	this->BtorGOverR0_f = new real_t[GetNr()+1];
    this->BpolIOverR0 = new real_t[GetNr()];
    this->BpolIOverR0_f = new real_t[GetNr()+1];
    this->rotTransf = new real_t[GetNr()];
    this->rotTransf_f = new real_t[GetNr()+1];
    //this->BflxKOverR0 = new real_t[GetNr()*ntheta_interp*nphi_interp]; // TODO: Is this needed in this form?
    //this->BflxKOverR0_f = new real_t[(GetNr()+1)(ntheta_interp)*(nphi_interp)]; // TODO: Is it needed in this form?

    this->acc_r = gsl_interp_accel_alloc();
    this->acc_theta = gsl_interp_accel_alloc();

    LoadMagneticFieldData(mf, frmt);
}

/**
 * Destructor.
 */
NumericBRadialGridGenerator::~NumericBRadialGridGenerator() {
    gsl_multimin_fminimizer_free(gsl_multi_fmin);

    if (this->rf_provided != nullptr)
        delete [] this->rf_provided;
	
	if (this->input_r != nullptr)
		delete [] this->input_r;
	
	if (this->psi != nullptr)
		delete [] this->psi;
	if (this->theta != nullptr)
		delete [] this->theta;
	if (this->phi != nullptr)
		delete [] this->phi;
	if (this->dataG != nullptr)
		delete [] this->dataG;
	if (this->dataI != nullptr)
		delete [] this->dataI;
	if (this->dataiota != nullptr)
		delete [] this->dataiota;
	//if (this->dataK != nullptr)
	//	delete [] this->dataK;
	if (this->dataBdotGradphi != nullptr)
		delete [] this->dataBdotGradphi;
	if (this->datagtt != nullptr)
		delete [] this->datagtt;
	if (this->datagtp!= nullptr)
		delete [] this->datagtp;
	if (this->dataB != nullptr)
		delete [] this->dataB;
	if (this->dataJacobian != nullptr)
		delete [] this->dataJacobian;

    if (this->spline_R != nullptr) {
        gsl_spline_free(this->spline_psi);
        gsl_spline_free(this->spline_G);
        gsl_spline_free(this->spline_I);
        gsl_spline_free(this->spline_iota);
        //gsl_spline2d_free(this->spline_B);
        //gsl_spline2d_free(this->spline_Jacobian);
        
        //delete this->interp_K;
        delete this->interp_B;
        delete this->interp_Jacobian;
        delete this->interp_BdotGradphi;
        delete this->interp_gtt;
        delete this->interp_gtp;
    }

    gsl_interp_accel_free(this->acc_theta);
    gsl_interp_accel_free(this->acc_r);
}


/**
 * Make sure that the given poloidal angle is on the interval
 * [0, 2pi], and if not, shift so that it is.
 */
real_t NumericBRadialGridGenerator::_angleBounded(const real_t angle) const {
    if (angle >= 0 && angle <= 2*M_PI)
        return angle;

    // if angle > 0
    //   floor(...) = N
    // if angle < 0
    //   floor(...) = -(N+1)
    return angle - 2*M_PI*std::floor(angle/(2*M_PI));
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
        case FILE_FORMAT_DESC:
            d = LoadNumericBFromDESC(sf);
            break;

        default:
            throw FVMException(
                "NumericBRadialGrid: Unrecognized file format specified for magnetic field data: %d.",
                frmt
            );
    }

    auto convert_data = [](const double *a, const len_t nx, const len_t ny, const len_t nz) {
        real_t *d = new real_t[nx*ny*nz];
        for (len_t i = 0; i < nx*ny*nz; i++)
            d[i] = (real_t)a[i];

        return d;
    };

	this->npsi = d->npsi;
	this->ntheta = d->ntheta;
	this->nphi = d->nphi;

	this->R0 = d->R0;

    // Set magnetic field data
    this->psi      = convert_data(d->psi, npsi, 1, 1);
    this->theta    = convert_data(d->theta, ntheta, 1, 1);
    this->phi      = convert_data(d->phi, nphi, 1, 1);
    
    this->dataG    = convert_data(d->G, npsi, 1, 1);
    this->dataI    = convert_data(d->I, npsi, 1, 1); 

    //this->dataK    = convert_data(d->K, npsi, ntheta, nphi);
    this->dataBdotGradphi = convert_data(d->BdotGradphi, npsi, ntheta, nphi);

    this->dataB        = convert_data(d->B, npsi, ntheta, 1); 
    this->dataJacobian = convert_data(d->Jacobian, npsi, ntheta, 1);

    this->datagtt      = convert_data(d->gtt, npsi, ntheta, nphi); // TODO: Correct dim?
    this->datagtp      = convert_data(d->gtp, npsi, ntheta, nphi); // TODO: Correct dim?

	// TODO: extrapolate to r=0 and r=a?

	// Ensure that theta=0 exists...
	if (this->theta[0] != 0)
		throw FVMException("NumericBRadialGrid: All numerical data must be given in theta = 0.");
    

    // TODO: Update this depending on how data looks!
	// Add point at theta=2*pi if it does not already exist...
	if (this->theta[this->ntheta-1] != 2*M_PI) {
		this->theta = this->addThetaDataPoint(this->theta, 1, this->ntheta, 1);
		this->theta[this->ntheta] = 2*M_PI;
		//this->dataK = this->addThetaDataPoint(this->dataK, this->npsi, this->ntheta, this->nphi);
		this->dataBdotGradphi = this->addThetaDataPoint(this->dataBdotGradphi, this->npsi, this->ntheta, this->nphi);
		this->dataB = this->addThetaDataPoint(this->dataB, this->npsi, this->ntheta, 1);
        this->dataJacobian = this->addThetaDataPoint(this->dataJacobian, this->npsi, this->ntheta, 1);
        this->datagtt = this->addThetaDataPoint(this->datagtt, this->npsi, this->ntheta, this->nphi);
        this->datagtp = this->addThetaDataPoint(this->datagtp, this->npsi, this->ntheta, this->nphi);

		// r grid has now been extended...
		this->ntheta++;
	} else {
		// Ensure that values at theta=0 and theta=2*pi are identical
		for (len_t ir = 0; ir < this->npsi; ir++) {
			len_t idx0 = ir*this->ntheta + 0;
			len_t idx1 = ir*this->ntheta + (this->ntheta-1);

            if (this->dataB[idx0] != this->dataB[idx1])
				this->dataB[idx1] = this->dataB[idx0];
			if (this->dataJacobian[idx0] != this->dataJacobian[idx1])
				this->dataJacobian[idx1] = this->dataJacobian[idx0];

            for (len_t iphi; iphi < this->nphi; iphi++){
                idx0 = (ir*this->ntheta + 0)*this->nphi + this->iphi;
                idx1 = (ir*this->ntheta + (this->ntheta-1))*this->nphi + this->iphi;
                
                //if (this->dataK[idx0] != this->dataK[idx1])
                //    this->dataK[idx1] = this->dataK[idx0];
                if (this->datagtt[idx0] != this->datagtt[idx1])
                    this->datagtt[idx1] = this->datagtt[idx0];
                if (this->datagtp[idx0] != this->datagtp[idx1])
                    this->datagtp[idx1] = this->datagtp[idx0];

            }

		}
	}

    // Evaluate minor radius in outer midplane
    this->input_r = new real_t[this->npsi];
    for (len_t i = 0; i < this->npsi; i++)
        this->input_r[i] = 0; // TODO!

    // Verify that the magnetic field contains exactly one maximum and
    // one minimum

    // TODO: Should it only contain one of each?
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
 * Extend the given array with one element at theta=2*pi.
 */
real_t *NumericBRadialGridGenerator::addThetaDataPoint(
	const real_t *x, const len_t npsi, const len_t ntheta, const len_t nphi
) {
	real_t *arr = new real_t[npsi*(ntheta+1)*nphi];

	for (len_t ir = 0, i, j, j_temp; ir < npsi; ir++){
        for (len_t ithetaphi = 0; i < nphi; i++, j++){
            arr[i] = x[j];
            arr[i+ntheta*nphi] = x[j];
        }
        for (len_t ithetaphi = nphi; i < ntheta*nphi; i++, j++){
            arr[i] = x[j];
        }
        i += nphi;
    }
    
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
    this->spline_G    = gsl_spline_alloc(gsl_interp_steffen, this->npsi); 
    this->spline_I    = gsl_spline_alloc(gsl_interp_steffen, this->npsi);

    const gsl_interp2d_type *splineType = gsl_interp2d_bicubic; //or ..._bilinear
    //this->spline_B           = gsl_spline2d_alloc(splineType, this->npsi, this->ntheta);
    //this->spline_Jacobian    = gsl_spline2d_alloc(splineType, this->npsi, this->ntheta);

    gsl_spline_init(this->spline_psi, this->input_r, this->psi,   this->npsi);
    gsl_spline_init(this->spline_G,   this->input_r, this->dataG, this->npsi);
    gsl_spline_init(this->spline_I,   this->input_r, this->dataI, this->npsi); 

    //gsl_spline2d_init(this->spline_B,        this->input_r, this->theta, this->dataB, this->npsi, this->ntheta);
    //gsl_spline2d_init(this->spline_Jacobian, this->input_r, this->theta, this->dataJacobian, this->npsi, this->ntheta);
    
    enum FVM::Interpolator3D::interp_method interp_meth = FVM::Interpolator3D::INTERP_LINEAR; // TODO: possibly implement cubic?
    //this->interp_K           = new FVM::Interpolator3D(this->npsi, this->ntheta, this->nphi, this->input_r, this->theta, this->phi, this->dataK, nullptr, interp_meth);
    this->interp_B           = new FVM::Interpolator3D(this->npsi, this->ntheta, this->nphi, this->input_r, this->theta, this->phi, this->dataB, nullptr, interp_meth);
    this->interp_Jacobian    = new FVM::Interpolator3D(this->npsi, this->ntheta, this->nphi, this->input_r, this->theta, this->phi, this->dataJacobian, nullptr, interp_meth);
    this->interp_BdotGradphi = new FVM::Interpolator3D(this->npsi, this->ntheta, this->nphi, this->input_r, this->theta, this->phi, this->dataBdotGradphi, nullptr, interp_meth);
    this->interp_gtt         = new FVM::Interpolator3D(this->npsi, this->ntheta, this->nphi, this->input_r, this->theta, this->phi, this->datagtt, nullptr, interp_meth);
    this->interp_gtp         = new FVM::Interpolator3D(this->npsi, this->ntheta, this->nphi, this->input_r, this->theta, this->phi, this->datagtp, nullptr, interp_meth);

	// Reference quantities
	for (len_t i = 0; i < GetNr(); i++) {
		real_t
			BtorG = gsl_spline_eval(
					this->spline_G, this->r[i], this->acc_r
				),
			BpolI = gsl_spline_eval(
					this->spline_I, this->r[i], this->acc_r
				), 
			iota = gsl_spline_eval(
					this->spline_iota, this->r[i], this->acc_r
				), 
            psip = gsl_spline_eval_deriv(
                    this->spline_psi, this->r[i], this->acc_r
                );

		this->BtorGOverR0[i] = BtorG / this->R0;
		this->BpolIOverR0[i] = BpolI / this->R0; 
		this->rotTransf[i]   = iota; 
        this->psiPrimeRef[i] = psip; // TODO: We could use the toroidal current here, and be consistent with stellarators. Would mean using d psi_t/ d r in bootstrap
	}
	for (len_t i = 0; i < GetNr()+1; i++) {
		real_t
			BtorG = gsl_spline_eval(
					this->spline_G, this->r_f[i], this->acc_r
				),
			BpolI = gsl_spline_eval(
					this->spline_I, this->r_f[i], this->acc_r
				), 
			iota = gsl_spline_eval(
					this->spline_iota, this->r_f[i], this->acc_r
				), 
            psip = gsl_spline_eval_deriv(
                    this->spline_psi, this->r_f[i], this->acc_r
                );

		this->BtorGOverR0_f[i] = BtorG / this->R0;
		this->BpolIOverR0_f[i] = BpolI / this->R0; 
		this->rotTransf_f[i]   = iota; 
        this->psiPrimeRef_f[i] = psip; // TODO: Does this even make sense?
	}
	
	rGrid->SetReferenceMagneticFieldData(
		this->BtorGOverR0, this->BtorGOverR0_f,
		this->BpolIOverR0, this->BpolIOverR0_f,
		this->psiPrimeRef, this->psiPrimeRef_f,
		this->R0
	);
    SetStellaratorData(this->rotTransf, this->rotTransf_f);

    this->isBuilt = true;
    return true;
}

/**
 * Calculate the Jacobian J, normalized to R0, at the specified radius and
 * poloidal angle.
 *
 * r:     Minor radius.
 * theta: Poloidal angle.
 */
real_t NumericBRadialGridGenerator::JacobianAtThetaPhi(
    const real_t radius, const real_t theta, const real_t phi
) {
    
    real_t t[1] = {this->_angleBounded(theta)};
    real_t p[1] = {this->_angleBounded(phi)};
    real_t r[1] = {radius}

    real_t *Jacobian = new real_t[1];

    this->interp_Jacobian->Eval(1,1,1, &r, &t, &p, nullptr, Jacobian);

    return Jacobian[0] / this->R0;
}

/**
 * Calculate B\dot Grad phi at the given poloidal and toroidal angle.
 *
 * r:     Minor radius.
 * theta: Poloidal angle.
 * phi:   Toroidal angle
 */
real_t NumericBRadialGridGenerator::BdotGradphiAtThetaPhi(
    const real_t radius, const real_t theta, const real_t phi
) {
    real_t t[1] = {this->_angleBounded(theta)};
    real_t p[1] = {this->_angleBounded(phi)};
    real_t r[1] = {radius}

    real_t *BdotGradphi = new real_t[1];

    this->interp_BdotGradphi->Eval(1,1,1, &r, &t, &p, nullptr, BdotGradphi);

	return BdotGradphi[0];
}

/**
 * Calculate g_{\theta\theta}/J^2 at the given poloidal and toroidal angle.
 *
 * r:     Minor radius.
 * theta: Poloidal angle.
 * phi:   Toroidal angle
 */
real_t NumericBRadialGridGenerator::gttAtThetaPhi(
    const real_t radius, const real_t theta, const real_t phi
) {
    real_t t[1] = {this->_angleBounded(theta)};
    real_t p[1] = {this->_angleBounded(phi)};
    real_t r[1] = {radius}

    real_t *gtt = new real_t[1];
    
    this->interp_gtt->Eval(1,1,1, &r, &t, &p, nullptr, gtt);

    real_t J = JacobianAtThetaPhi(radius, theta, phi);
	return gtt[0] / (J*J) / (R0 * R0);
}

/**
 * Calculate g_{\theta\varphi}/J^2 at the given poloidal and toroidal angle.
 *
 * r:     Minor radius.
 * theta: Poloidal angle.
 * phi:   Toroidal angle
 */
real_t NumericBRadialGridGenerator::gtpAtThetaPhi(
    const real_t radius, const real_t theta, const real_t phi
) {
    real_t t[1] = {this->_angleBounded(theta)};
    real_t p[1] = {this->_angleBounded(phi)};
    real_t r[1] = {radius}

    real_t *gtp = new real_t[1];

    this->interp_gtp->Eval(1,1,1, &r, &t, &p, nullptr, gtp);

    real_t J = JacobianAtThetaPhi(radius, theta, phi);
	return gtp[0] / (J*J) / (R0 * R0);
}

/**
 * Evaluate all the geometric quantities in one go.
 */
void NumericBRadialGridGenerator::EvaluateGeometricQuantities(
    const real_t r, const real_t theta, real_t phi, real_t &B, real_t &Jacobian, real_t &BdotGradphi, real_t &gttOverJ2, real_t &gtpOverJ2
) {
    real_t t = this->_angleBounded(theta);
    
    real_t p = this->_angleBounded(phi);

    Jacobian = JacobianAtThetaPhi(r, t, p);
    
    B = BAtThetaPhi(r, t, p);

    BdotGradphi = BdotGradphiAtThetaPhi(r, t, p);

    gttOverJ2 = gttAtThetaPhi(r, t, p);

    gtpOverJ2 = gtpAtThetaPhi(r, t, p);
}

/**
 * Evaluate magnetic field strength B at given poloidal angle
 * and radius.
 */
real_t NumericBRadialGridGenerator::EvalB(const real_t radius, const real_t theta, const real_t phi) {
    
    real_t t[1] = {this->_angleBounded(theta)};
    real_t p[1] = {this->_angleBounded(phi)};
    real_t r[1] = {radius}

    real_t *B = new real_t[1];

    this->interp_B->Eval(1,1,1, &r, &t, &p, nullptr, B);
    
    return B[0];
}

real_t NumericBRadialGridGenerator::BAtThetaPhi(const len_t ir, const real_t theta, const real_t phi) {
	return EvalB(this->r[ir], theta, phi);
}
real_t NumericBRadialGridGenerator::BAtThetaPhi_f(const len_t ir, const real_t theta, const real_t phi) {
	return EvalB(this->r_f[ir], theta, phi);
}

// The remaining functions are related to determining theta_Bmin and theta_Bmax
// with a gsl fmin algorithm
struct EvalBParams {len_t ir; NumericBRadialGridGenerator* rgg; int_t sgn;};
real_t gslEvalB(const gsl_vector *v, void *par){
    real_t theta = gsl_vector_get(v, 0);
    real_t phi = gsl_vector_get(v, 1);
    EvalBParams *params = (EvalBParams *) par;
    return params->sgn*params->rgg->BAtThetaPhi(params->ir,theta,phi);
}
real_t gslEvalB_f(real_t theta, void *par){
    EvalBParams *params = (EvalBParams *) par;
    return params->sgn*params->rgg->BAtThetaPhi_f(params->ir,theta, phi);
}


// Tolerances and max number of iterations for gsl magnetic field minimizer
const len_t MAX_NUM_ITER = 30;
const real_t EPSABS = 1e-6;
const real_t EPSREL = 0;
const real_t STEP = 2*M_PI / 100; // TODO: Ok? How many wiggles can we expect at the most?
/**
 * Finds the extremum of the magnetic field on the interval [0,2*pi]. 
 * If sgn=1, returns the minimum.
 * If sgn=-1, returns the maximum.
 */
real_t NumericBRadialGridGenerator::FindMagneticFieldExtremum(
    len_t ir, int_t sgn, fluxGridType fluxGridType
) {
    real_t theta_lim_lower = 0, theta_lim_upper = 2*M_PI;
    
    real_t theta_guess = 0, phi_guess = 0;
    real_t B_opt = sgn*std::numeric_limits<real_t>::infinity();
    real_t B;

    for (real_t theta=0; theta<2*M_PI; theta+=STEP){
        for (real_t phi=0; phi<2*M_PI; phi+=STEP){
            if (fluxGridType == FLUXGRIDTYPE_DISTRIBUTION) {
                B = BAtThetaPhi(ir, theta, phi);
            } else {
                B = BAtThetaPhi_f(ir, theta, phi);
            }
            if (sgn*B < sgn*B_opt){
                B_opt = B;
                theta_guess = theta;
                phi_guess = phi;
            }
        }
    }
    '''
    real_t theta_lim_lower = std::max(0.0, theta_guess - STEP);
    real_t theta_lim_upper = std::min(2*M_PI, theta_guess + STEP);
    real_t phi_lim_lower = std::max(0.0, phi_guess - STEP);
    real_t phi_lim_upper = std::min(2*M_PI, phi_guess + STEP);
    '''
    gsl_vector *guess = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, theta_guess);
    gsl_vector_set(x, 1, phi_guess);

    gsl_vector *step = gsl_vector_alloc(2);
    gsl_vector_set_all(ss, STEP / 2); // TODO: OK Step size?
	
    EvalBParams params = {ir, this, sgn};
    gsl_multimin_function gsl_func;
    gsl_func.n = 2;
    if(fluxGridType == FLUXGRIDTYPE_DISTRIBUTION)
        gsl_func.f = &(gslEvalB); 
    else
        gsl_func.f = &(gslEvalB_f); 
    gsl_func.params = &(params); // TODO: Should this be without &?

    
    // otherwise, find extremum with fmin algorithm
    gsl_multimin_fminimizer_set( 
        gsl_multi_fmin, &gsl_func, guess, step, 
    );

    int status;
    for(len_t iter=0; iter<MAX_NUM_ITER; iter++){
        gsl_multimin_fminimizer_iterate(gsl_multi_fmin);

        size = gsl_multimin_fminimizer_size(gsl_multi_fmin);
        status = gsl_multimin_test_size(size, EPSABS);
        if(status == GSL_SUCCESS)
            break;
    }

    gsl_vector_free(guess);
    gsl_vector_free(step);

    real_t theta = gsl_vector_get(gsl_multi_fmin->x, 0);
    real_t phi   = gsl_vector_get(gsl_multi_fmin->x, 1);
  
    real_t extremum = theta; // TODO: We disregard phi, right?
    if(extremum < 2*EPSABS || extremum > 2*M_PI - 2*EPSABS)
        return 0;
    else if (fabs(M_PI-extremum) < 2*EPSABS)
        return M_PI;
    else
        return extremum; 
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

/**
 * Returns a list of toroidal angles on which the flux
 * surfaces returned by 'GetFluxSurfaceR()' and 'GetFluxSurfaceZ()'
 * are defined.
 */
const real_t *NumericBRadialGridGenerator::GetToroidalAngle() {
	real_t *phi = new real_t[nphi];
	for (len_t i = 0; i < nphi; i++)
		phi[i] = this->phi[i];
	
	return phi;
}


/** TODO: Is this good?
 * Save the magnitude of the magnetic field vector to the named
 * output file (saved using the 'SFile' API).
 *
 * filename: Name of file to save data to.
 */
void NumericBRadialGridGenerator::__SaveB(const char *filename) {
    // DEBUG: Save magnetic field
    const len_t NTHETA = 100, NPHI = 100;
    real_t **B = new real_t*[GetNr()];
    B[0] = new real_t[GetNr()*NTHETA];
    for (len_t i = 0; i < GetNr(); i++) {
        if (i > 0)
            B[i] = B[i-1] + NTHETA;

        for (len_t j = 0; j < NPHI; j++)
            for (len_t k = 0; k < NTHETA; k++)
                B[i][j*NTHETA + k] = this->BAtThetaPhi(i, k*2*M_PI/NTHETA, j*2*M_PI/NPHI);
    }

    SFile *sf = SFile::Create(filename, SFILE_MODE_WRITE);
    sf->WriteArray("B", B, GetNr(), NTHETA);
    sf->Close();
}

