/**
 * Implementation of the numeric magnetic field radial grid generator. This
 * grid generator loads a numeric magnetic field from the specified file and
 * builds a correspondingly shaped radial grid.
 */
#include <algorithm>
#include <gsl/gsl_interp.h>
#include "FVM/Grid/Stellarator/NumericStellaratorRadialGridGenerator.hpp"


using namespace DREAM::FVM;


/**
 * Constructor for a uniform radial grid.
 *
 * nr:            Number of radial grid points in uniform radial *distribution* grid. 
 * r0:            Value of innermost point on radial *flux* grid.
 * ra:            Value of outermost point on radial *flux* grid.
 * ntheta_interp: Poloidal angle resolution in quadrature for flux surface and bounce averages.
 * nphi_interp:   Toroidal angle resolution in quadrature for flux surface and bounce averages.
 */
NumericStellaratorRadialGridGenerator::NumericStellaratorRadialGridGenerator(
    const len_t nr, const real_t r0, const real_t ra,
    const real_t R0, const len_t nfp, struct eq_data* eqdata,
	const len_t ntheta_interp, const len_t nphi_interp
) : RadialGridGeneratorStellarator(nr), providedData(eqdata), rMin(r0), rMax(ra) {

    this->nfp = nfp;
    this->R0 = R0;

    this->nrho = eqdata->nrho;
    this->ntheta = eqdata->ntheta;
    this->nphi = eqdata->nphi;
	this->ntheta_interp = ntheta_interp;
	this->nphi_interp   = nphi_interp;

	this->BtorGOverR0 = new real_t[GetNr()];
	this->BtorGOverR0_f = new real_t[GetNr()+1];
    this->BpolIOverR0 = new real_t[GetNr()];
    this->BpolIOverR0_f = new real_t[GetNr()+1];
    this->iota = new real_t[GetNr()];
    this->iota_f = new real_t[GetNr()+1];
    //this->BflxKOverR0 = new real_t[GetNr()*ntheta_interp*nphi_interp]; 
    //this->BflxKOverR0_f = new real_t[(GetNr()+1)(ntheta_interp)*(nphi_interp)]; 

    this->acc_r = gsl_interp_accel_alloc();
    this->acc_theta = gsl_interp_accel_alloc();
}

/**
 * Constructor for a custom radial grid.
 *
 * r_f:           Radial flux grid desired.
 * nr:            Number of points on distribution grids corresponding to 
 *                'r_f' (i.e. the number of points in 'r_f', minus one).
 * ntheta_interp: Poloidal angle resolution in quadrature for flux surface and bounce averages.
 * nphi_interp:   Toroidal angle resolution in quadrature for flux surface and bounce averages.
 */
NumericStellaratorRadialGridGenerator::NumericStellaratorRadialGridGenerator(
    const real_t *r_f, const len_t nr,
    const real_t R0, const len_t nfp, struct eq_data* eqdata,
	const len_t ntheta_interp, const len_t nphi_interp
) : RadialGridGeneratorStellarator(nr), providedData(eqdata) {

    this->nfp = nfp;
    this->R0 = R0;

    this->rf_provided = new real_t[nr+1];
    for (len_t i = 0; i < nr+1; i++)
        this->rf_provided[i] = r_f[i];
    
    this->nrho = eqdata->nrho;
    this->ntheta = eqdata->ntheta;
    this->nphi = eqdata->nphi;

	this->ntheta_interp = ntheta_interp;
	this->nphi_interp   = nphi_interp;

	this->BtorGOverR0 = new real_t[GetNr()];
	this->BtorGOverR0_f = new real_t[GetNr()+1];
    this->BpolIOverR0 = new real_t[GetNr()];
    this->BpolIOverR0_f = new real_t[GetNr()+1];
    this->iota = new real_t[GetNr()];
    this->iota_f = new real_t[GetNr()+1];
    //this->BflxKOverR0 = new real_t[GetNr()*ntheta_interp*nphi_interp]; 
    //this->BflxKOverR0_f = new real_t[(GetNr()+1)(ntheta_interp)*(nphi_interp)]; 

    this->acc_r = gsl_interp_accel_alloc();
    this->acc_theta = gsl_interp_accel_alloc();
}

/**
 * Destructor.
 */
NumericStellaratorRadialGridGenerator::~NumericStellaratorRadialGridGenerator() {

    if (this->rf_provided != nullptr)
        delete [] this->rf_provided;
	

    if (this->spline_G != nullptr) {
        gsl_spline_free(this->spline_G);
        gsl_spline_free(this->spline_I);
        gsl_spline_free(this->spline_iota);
        gsl_spline_free(this->spline_psi);
        
        //delete this->interp_K;
        delete this->interp_B;
        delete this->interp_Jacobian;
        delete this->interp_BdotGradphi;
        delete this->interp_gtt;
        delete this->interp_gtp;
        delete this->interp_lambdat;
        delete this->interp_lambdap;
    }

    gsl_interp_accel_free(this->acc_theta);
    gsl_interp_accel_free(this->acc_r);

	delete this->providedData; // Is this ok?
}


/**
 * Make sure that the given poloidal angle is on the interval
 * [0, 2pi], and if not, shift so that it is.
 */
real_t NumericStellaratorRadialGridGenerator::_angleBounded(const real_t angle) const {
    if (angle >= 0 && angle <= 2*M_PI)
        return angle;

    // if angle > 0
    //   floor(...) = N
    // if angle < 0
    //   floor(...) = -(N+1)
    return angle - 2*M_PI*std::floor(angle/(2*M_PI));
}

/**
 * Rebuild this numeric radial B grid.
 */
bool NumericStellaratorRadialGridGenerator::Rebuild(const real_t, RadialGridStellarator *rGrid) {
    this->r   = new real_t[GetNr()];
    this->r_f = new real_t[GetNr()+1];

    real_t
        *dr   = new real_t[GetNr()],
        *dr_f = new real_t[GetNr()-1];

    // If a custom grid has been specified, set it here,
    // otherwise we generate a uniform radial grid
    if (rf_provided == nullptr) {
        if (rMin < 0)
            throw FVMException("StellaratorRadialGrid: rMin < 0.");
        else if (rMax > this->providedData->rho[this->nrho-1])
            throw FVMException(
                "StellaratorRadialGrid: Maximum r available in numeric magnetic "
                "field data is rMax = %.3f, but r = %.3f is required for radial "
                "grid. (Delta: %e)", this->providedData->rho[this->nrho-1], rMax,
                std::abs(this->providedData->rho[this->nrho-1]/rMax-1)
            );
        else if (rMin >= rMax)
            throw FVMException("StellaratorRadialGrid: rMin must be strictly less than rMax.");

        for (len_t i = 0; i < GetNr(); i++)
            dr[i] = (rMax-rMin) / GetNr();
        for (len_t i = 0; i < GetNr()+1; i++)
            r_f[i] = rMin + i*dr[0];
    } else {
        if (rf_provided[0] < 0)
            throw FVMException("StellaratorRadialGrid: First point on custom radial grid is less than zero.");
        else if (rf_provided[GetNr()] > this->providedData->rho[this->nrho-1])
            throw FVMException(
                "StellaratorRadialGrid: Last point on custom radial grid may not be greater than "
                "the maximum r available in the numeric magnetic field data, rMax = %.3f.",
                this->providedData->rho[this->nrho-1]
            );
        else if (rf_provided[0] >= rf_provided[GetNr()])
            throw FVMException("StellaratorRadialGrid: The first point on the custom radial grid must be strictly less than the last point.");

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
    this->spline_psi  = gsl_spline_alloc(gsl_interp_steffen, this->nrho);
    this->spline_G    = gsl_spline_alloc(gsl_interp_steffen, this->nrho); 
    this->spline_I    = gsl_spline_alloc(gsl_interp_steffen, this->nrho);
    this->spline_iota = gsl_spline_alloc(gsl_interp_steffen, this->nrho);

    gsl_spline_init(this->spline_G,    this->providedData->rho, this->providedData->dataG, this->nrho);
    gsl_spline_init(this->spline_I,    this->providedData->rho, this->providedData->dataI, this->nrho); 
    gsl_spline_init(this->spline_iota, this->providedData->rho, this->providedData->dataiota, this->nrho);
    gsl_spline_init(this->spline_psi,  this->providedData->rho, this->providedData->datapsi, this->nrho);
    
    enum FVM::Interpolator3DSpatial::interp_method interp_meth = FVM::Interpolator3DSpatial::INTERP_LINEAR; 
    //this->interp_K           = new FVM::Interpolator3DSpatial(this->nphi, this->nrho, this->ntheta, this->providedData->phi, this->providedData->rho, this->providedData->theta, this->providedData->dataK,interp_meth);
    this->interp_B           = new FVM::Interpolator3DSpatial(this->nphi, this->nrho, this->ntheta, this->providedData->phi, this->providedData->rho, this->providedData->theta, this->providedData->dataB, interp_meth);
    this->interp_Jacobian    = new FVM::Interpolator3DSpatial(this->nphi, this->nrho, this->ntheta, this->providedData->phi, this->providedData->rho, this->providedData->theta, this->providedData->dataJacobian,interp_meth);
    this->interp_BdotGradphi = new FVM::Interpolator3DSpatial(this->nphi, this->nrho, this->ntheta, this->providedData->phi, this->providedData->rho, this->providedData->theta, this->providedData->dataBdotGradphi,interp_meth);
    this->interp_gtt         = new FVM::Interpolator3DSpatial(this->nphi, this->nrho, this->ntheta, this->providedData->phi, this->providedData->rho, this->providedData->theta, this->providedData->datagtt,interp_meth);
    this->interp_gtp         = new FVM::Interpolator3DSpatial(this->nphi, this->nrho, this->ntheta, this->providedData->phi, this->providedData->rho, this->providedData->theta, this->providedData->datagtp,interp_meth);
    this->interp_lambdat     = new FVM::Interpolator3DSpatial(this->nphi, this->nrho, this->ntheta, this->providedData->phi, this->providedData->rho, this->providedData->theta, this->providedData->datalambdat,interp_meth);
    this->interp_lambdap     = new FVM::Interpolator3DSpatial(this->nphi, this->nrho, this->ntheta, this->providedData->phi, this->providedData->rho, this->providedData->theta, this->providedData->datalambdap,interp_meth);

	// Reference quantities
	for (len_t i = 0; i < GetNr(); i++) {
		real_t
			BtorG = gsl_spline_eval(
					this->spline_G, this->r[i], this->acc_r
				),
			BpolI = gsl_spline_eval(
					this->spline_I, this->r[i], this->acc_r
				), 
			rottr = gsl_spline_eval(
					this->spline_iota, this->r[i], this->acc_r
				), 
            psip = gsl_spline_eval_deriv(
                    this->spline_psi, this->r[i], this->acc_r
                );

		this->BtorGOverR0[i] = BtorG / this->R0;
		this->BpolIOverR0[i] = BpolI / this->R0; 
		this->iota[i]        = rottr; 
        this->psiPrimeRef[i] = psip; 
	}
	for (len_t i = 0; i < GetNr()+1; i++) {
		real_t
			BtorG = gsl_spline_eval(
					this->spline_G, this->r_f[i], this->acc_r
				),
			BpolI = gsl_spline_eval(
					this->spline_I, this->r_f[i], this->acc_r
				), 
			rottr = gsl_spline_eval(
					this->spline_iota, this->r_f[i], this->acc_r
				), 
            psip = gsl_spline_eval_deriv(
                    this->spline_psi, this->r_f[i], this->acc_r
                );

		this->BtorGOverR0_f[i] = BtorG / this->R0;
		this->BpolIOverR0_f[i] = BpolI / this->R0; 
		this->iota_f[i]   = rottr; 
        this->psiPrimeRef_f[i] = psip; 
	}
	
	rGrid->SetReferenceMagneticFieldData(
		this->BtorGOverR0, this->BtorGOverR0_f,
		this->BpolIOverR0, this->BpolIOverR0_f,
		this->psiPrimeRef, this->psiPrimeRef_f,
		this->R0
	);
    rGrid->SetStellaratorData(this->iota, this->iota_f);

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
real_t NumericStellaratorRadialGridGenerator::JacobianAtThetaPhi(
    const real_t radius, const real_t theta, const real_t phi
) {
    
    real_t t[1] = {this->_angleBounded(theta)};
    real_t p[1] = {this->_angleBounded(phi)};
    real_t r[1] = {radius};

    real_t *Jacobian = new real_t[1];

    this->interp_Jacobian->Eval(1,1,1, p, r, t, Jacobian);

    return Jacobian[0] / this->R0;
}

/**
 * Calculate B\dot Grad phi at the given poloidal and toroidal angle.
 *
 * r:     Minor radius.
 * theta: Poloidal angle.
 * phi:   Toroidal angle
 */
real_t NumericStellaratorRadialGridGenerator::BdotGradphiAtThetaPhi(
    const real_t radius, const real_t theta, const real_t phi
) {
    real_t t[1] = {this->_angleBounded(theta)};
    real_t p[1] = {this->_angleBounded(phi)};
    real_t r[1] = {radius};

    real_t *BdotGradphi = new real_t[1];

    this->interp_BdotGradphi->Eval(1,1,1, p, r, t, BdotGradphi);

	return BdotGradphi[0] * this->R0;
}

/**
 * Calculate g_{\theta\theta}/J^2 at the given poloidal and toroidal angle.
 *
 * r:     Minor radius.
 * theta: Poloidal angle.
 * phi:   Toroidal angle
 */
real_t NumericStellaratorRadialGridGenerator::gttAtThetaPhi(
    const real_t radius, const real_t theta, const real_t phi
) {
    real_t t[1] = {this->_angleBounded(theta)};
    real_t p[1] = {this->_angleBounded(phi)};
    real_t r[1] = {radius};

    real_t *gtt = new real_t[1];
    
    this->interp_gtt->Eval(1,1,1, p, r, t, gtt);

    real_t J = JacobianAtThetaPhi(radius, theta, phi);
	return gtt[0] / (J*J);
}

/**
 * Calculate g_{\theta\varphi}/J^2 at the given poloidal and toroidal angle.
 *
 * r:     Minor radius.
 * theta: Poloidal angle.
 * phi:   Toroidal angle
 */
real_t NumericStellaratorRadialGridGenerator::gtpAtThetaPhi(
    const real_t radius, const real_t theta, const real_t phi
) {
    real_t t[1] = {this->_angleBounded(theta)};
    real_t p[1] = {this->_angleBounded(phi)};
    real_t r[1] = {radius};

    real_t *gtp = new real_t[1];
    real_t *gtt = new real_t[1];
    real_t *lambdat = new real_t[1];
    real_t *lambdap = new real_t[1];
    
    this->interp_gtp->Eval(1,1,1, p, r, t, gtp);
    this->interp_gtt->Eval(1,1,1, p, r, t, gtt);
    this->interp_lambdat->Eval(1,1,1, p, r, t, lambdat);
    this->interp_lambdap->Eval(1,1,1, p, r, t, lambdap);

    real_t J = JacobianAtThetaPhi(radius, theta, phi);
	return (gtp[0] * (1 + lambdat[0]) - gtt[0] * lambdap[0]) / (J*J) ;
}

/**
 * Evaluate all the geometric quantities in one go.
 */
void NumericStellaratorRadialGridGenerator::EvaluateGeometricQuantities(
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
real_t NumericStellaratorRadialGridGenerator::EvalB(const real_t radius, const real_t theta, const real_t phi) {
    
    real_t t[1] = {this->_angleBounded(theta)};
    real_t p[1] = {this->_angleBounded(phi)};
    real_t r[1] = {radius};

    real_t *B = new real_t[1];

    this->interp_B->Eval(1,1,1, p, r, t, B);
    
    return B[0];
}

real_t NumericStellaratorRadialGridGenerator::BAtThetaPhi(const len_t ir, const real_t theta, const real_t phi) {
	return EvalB(this->r[ir], theta, phi);
}
real_t NumericStellaratorRadialGridGenerator::BAtThetaPhi_f(const len_t ir, const real_t theta, const real_t phi) {
	return EvalB(this->r_f[ir], theta, phi);
}


/**
 * Returns a list of poloidal angles on which the flux
 * surfaces returned by 'GetFluxSurfaceR()' and 'GetFluxSurfaceZ()'
 * are defined.
 */
const real_t *NumericStellaratorRadialGridGenerator::GetPoloidalAngle() {
	const len_t ntheta = this->GetNTheta();
	real_t *theta = new real_t[ntheta];

	for (len_t i = 0; i < ntheta; i++)
		theta[i] = 2*M_PI*i / ntheta;
	
	return theta;
}

/**
 * Returns a list of toroidal angles on which the flux
 * surfaces returned by 'GetFluxSurfaceR()' and 'GetFluxSurfaceZ()'
 * are defined.
 */
const real_t *NumericStellaratorRadialGridGenerator::GetToroidalAngle() {
	const len_t nphi = this->GetNPhi();
	real_t *phi = new real_t[nphi];

    real_t factor = 1;
    if (this->nfp > 0)
        factor /= 2*nfp;

	for (len_t i = 0; i < nphi; i++)
		phi[i] = 2*M_PI*i / nphi * factor; // TODO: Ok? 
	
	return phi;
}


/** TODO: Is this good?
 * Save the magnitude of the magnetic field vector to the named
 * output file (saved using the 'SFile' API).
 *
 * filename: Name of file to save data to.
 */
void NumericStellaratorRadialGridGenerator::__SaveB(const char *filename) {
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

