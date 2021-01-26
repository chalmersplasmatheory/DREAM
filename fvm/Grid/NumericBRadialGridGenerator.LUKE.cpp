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
    d->npsi = fsize[1]==1 ? fsize[0] : fsize[1];

    // Poloidal angle coordinate grid
    d->theta = sf->GetList("equil/theta", fsize);
    d->ntheta = fsize[1]==1 ? fsize[0] : fsize[1];

    // Radial meshgrid
    d->R = sf->GetList("equil/ptx", fsize); ASSERT_DIMS("ptx");
    d->Z = sf->GetList("equil/pty", fsize); ASSERT_DIMS("pty");

    d->Br   = sf->GetList("equil/ptBx", fsize); ASSERT_DIMS("ptBx");
    d->Bz   = sf->GetList("equil/ptBy", fsize); ASSERT_DIMS("ptBz");
    d->Bphi = sf->GetList("equil/ptBPHI", fsize); ASSERT_DIMS("ptBPHI");

    return d;
}

