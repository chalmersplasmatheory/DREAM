/**
 * Routine for loading magnetic field data from a DESC
 * format data file.
 */

#include <string>
#include "FVM/Grid/NumericBRadialGridGenerator.hpp"

using namespace DREAM::FVM;


struct NumericBData *DREAM::FVM::LoadNumericBFromDESC(SFile *sf, const std::string&) {
    sfilesize_t fsize[3];

    struct NumericBData *d = new struct NumericBData;

    #define ASSERT_DIMS_1D(var) \
        if (fsize[0] != d->npsi)/* TODO: Correct?*/ \ 
            throw FVMException( \
                "%s: Invalid dimensions of vector '" var "' (%llu, %llu). " \
                "Expected (" LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT ").", \
                sf->filename.c_str(), fsize[0], fsize[1], d->ntheta, d->npsi \
            )
    #define ASSERT_DIMS_2D(var) \
        if (fsize[0] != d->npsi || fsize[1] != d->ntheta)/* TODO: Correct?*/ \ 
            throw FVMException( \
                "%s: Invalid dimensions of vector '" var "' (%llu, %llu). " \
                "Expected (" LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT ").", \
                sf->filename.c_str(), fsize[0], fsize[1], d->ntheta, d->npsi \
            )
    #define ASSERT_DIMS_3D(var) \
        if (fsize[0] != d->npsi || fsize[1] != d->ntheta || fsize[2] != d->nphi)/* TODO: Correct?*/ \ 
            throw FVMException( \
                "%s: Invalid dimensions of vector '" var "' (%llu, %llu). " \
                "Expected (" LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT ").", \
                sf->filename.c_str(), fsize[0], fsize[1], d->ntheta, d->npsi \
            )

    // Major radius // TODO: Change?
    d->R0 = (real_t)sf->GetScalar("equil/R0");

    // Poloidal flux coordinate grid
    d->psi = sf->GetList("equil/psi", fsize);
    //d->npsi = fsize[1]==1 ? fsize[0] : fsize[1];
	d->npsi = fsize[0]; // TODO: What does this do?

    // Poloidal angle coordinate grid
    d->theta = sf->GetList("equil/theta", fsize);
    //d->ntheta = fsize[1]==1 ? fsize[0] : fsize[1];
	d->ntheta = fsize[0]; // TODO: What does this do?

    // Toroidal angle coordinate grid
    d->phi = sf->GetList("equil/phi", fsize);
    //d->nphi = fsize[1]==1 ? fsize[0] : fsize[1];
	d->nphi = fsize[0]; // TODO: What does this do?

    // Radial meshgrid
    double **_G   = sf->GetDoubles("equil/G", fsize); ASSERT_DIMS_1D("G(psi)");
    double **_I   = sf->GetDoubles("equil/I", fsize); ASSERT_DIMS_1D("I(psi)");
    double **_iota = sf->GetDoubles("equil/iota", fsize); ASSERT_DIMS_1D("iota(psi)");
    //double **_K = sf->GetDoubles("equil/K", fsize); ASSERT_DIMS_3D("K(psi,theta,phi)");
    double **_BdotGradphi = sf->GetDoubles("equil/BdotGradphi", fsize); ASSERT_DIMS_3D("BdotGradphi(psi,theta,phi)");
    double **_B = sf->GetDoubles("equil/B", fsize); ASSERT_DIMS_2D("B(psi,theta)");
    double **_Jacobian = sf->GetDoubles("equil/Jacobian", fsize); ASSERT_DIMS_2D("Jacobian(psi,theta)");
    double **_gtt = sf->GetDoubles("equil/gtt", fsize); ASSERT_DIMS_3D("gtt(psi,theta,phi)");
    double **_gtp = sf->GetDoubles("equil/gtp", fsize); ASSERT_DIMS_3D("gtp(psi,theta,phi)");
    double **_lambdat = sf->GetDoubles("equil/lambdat", fsize); ASSERT_DIMS_3D("lambdat(psi,theta,phi)");
    double **_lambdap = sf->GetDoubles("equil/lambdap", fsize); ASSERT_DIMS_3D("lambdap(psi,theta,phi)");

	// We need only a pointer to the full chunk of
	// data, and so we free the 2D pointer here...
	d->G           = _G[0];             delete [] _G;
	d->I           = _I[0];             delete [] _I;
	d->iota        = _iota[0];          delete [] _iota;
	//d->K           = _K[0];             delete [] _K;
	d->BdotGradphi = _BdotGradphi[0];   delete [] _BdotGradphi;
	d->B           = _B[0];             delete [] _B;
	d->Jacobian    = _Jacobian[0];      delete [] _Jacobian;
    d->gtt         = _gtt[0];           delete [] _gtt;
	d->gtp         = _gtp[0];           delete [] _gtp;
	d->lambdat     = _lambdat[0];       delete [] _lambdat;
	d->lambdap     = _lambdap[0];       delete [] _lambdap;

    
    // TODO: Make DESC relevant
    // DREAM works with psi/R0 so we want to divide psi by R0 eventually
    for (len_t i = 0; i < d->npsi; i++)
        d->psi[i] /= d->R0;

    return d;
}

