/**
 * Routine for loading magnetic field data from a LUKE
 * format data file.
 */

#include <string>
#include "FVM/Grid/NumericBRadialGridGenerator.hpp"

using namespace DREAM::FVM;


struct NumericBData *DREAM::FVM::LoadNumericBFromLUKE(SFile *sf, const std::string&) {
    sfilesize_t fsize[2];

    struct NumericBData *d = new struct NumericBData;

    #define ASSERT_DIMS(var) \
        if (fsize[0] != d->ntheta || fsize[1] != d->npsi) \
            throw FVMException( \
                "%s: Invalid dimensions of vector '" var "' (%llu, %llu). " \
                "Expected (" LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT ").", \
                sf->filename.c_str(), fsize[0], fsize[1], d->ntheta, d->npsi \
            )

    d->name = sf->GetString("equil/id");

    // Magnetic axis coordinates
    d->Rp = (real_t)sf->GetScalar("equil/Rp");
    d->Zp = (real_t)sf->GetScalar("equil/Zp");

    // Poloidal flux coordinate grid
    d->psi = sf->GetList("equil/psi_apRp", fsize);
    //d->npsi = fsize[1]==1 ? fsize[0] : fsize[1];
	d->npsi = fsize[0];

    // Poloidal angle coordinate grid
    d->theta = sf->GetList("equil/theta", fsize);
    //d->ntheta = fsize[1]==1 ? fsize[0] : fsize[1];
	d->ntheta = fsize[0];

    // Radial meshgrid
    double **_R = sf->GetDoubles("equil/ptx", fsize); ASSERT_DIMS("ptx");
    double **_Z = sf->GetDoubles("equil/pty", fsize); ASSERT_DIMS("pty");

    double **_Br   = sf->GetDoubles("equil/ptBx", fsize); ASSERT_DIMS("ptBx");
    double **_Bz   = sf->GetDoubles("equil/ptBy", fsize); ASSERT_DIMS("ptBz");
    double **_Bphi = sf->GetDoubles("equil/ptBPHI", fsize); ASSERT_DIMS("ptBPHI");

	// We need only a pointer to the full chunk of
	// data, and so we free the 2D pointer here...
	d->R    = _R[0];    delete [] _R;
	d->Z    = _Z[0];    delete [] _Z;
	d->Br   = _Br[0];   delete [] _Br;
	d->Bz   = _Bz[0];   delete [] _Bz;
	d->Bphi = _Bphi[0]; delete [] _Bphi;

    return d;
}

