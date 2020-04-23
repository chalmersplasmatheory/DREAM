
#include "DREAM/IonHandler.hpp"

using namespace DREAM;


/**
 * Z is a list of atomic numbers of the ion species included (size nZ).
 * Example: a plasma containing deuterium and argon. 
 *          Z is size nZ=2 where Z[0] = 1, Z[1] = 18.
 * 
 * Densities are given on a nr x nzs grid, where nr is the number of radial grid points
 * and nzs the number of ion species (counting seperate charge states individually, 
 * where each ion has Z+1 possible charge states):
 *      nzs = sum_i (Z_i + 1), i = 0, 1, ..., nZ
 * The density of an ion of charge Z and charge state Z0 at radial point ir is given by
 *      n_i[nr * iz + ir]
 * and iz = offset(Z) + Z0, 
 */


IonHandler::IonHandler(FVM::RadialGrid *rg, FVM::UnknownQuantityHandler *u, len_t *Z, len_t NZ){
    rGrid = rg;
    unknowns = u;

    Zs  = Z;
    nr = rGrid->GetNr();
    nZ = NZ;

    niID = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

    Initialize();
}


IonHandler::~IonHandler(){
    DeallocateAll();
}

void IonHandler::Initialize(){
    nzs = 0;
    for (len_t it=0; it<nZ; it++){
        nzs += Zs[it]+1;
    }

    ZOffsets = new len_t[nZ];
    ZOffsets[0] = 0;
    if (nZ>0){
        for (len_t iz=1; iz<nZ; iz++){
            ZOffsets[iz] = ZOffsets[iz-1] + Zs[iz-1] + 1;
        }
    }
}


// Returns the density of ions which are characterised by 
// atomic number Z and charge number Z0 at radial index ir.
const real_t IonHandler::GetIonDensityAtZ(len_t ir, len_t Z, len_t Z0) const{
    real_t niReturn = 0;
    const real_t *n_i = unknowns->GetUnknownData(niID);
    for (len_t iz=0; iz<nZ; iz++)
        if (Zs[iz] == Z)
            niReturn += n_i[nr*(ZOffsets[iz]+Z0) + ir];
    

    return niReturn;
}

// Returns the density of ions which are characterised by 
// atomic number Z and charge number Z0 at radial index ir.
const real_t IonHandler::GetIonDensity(len_t ir, len_t iz, len_t Z0) const{
    const real_t *n_i = unknowns->GetUnknownData(niID);
    return n_i[nr*(ZOffsets[iz]+Z0) + ir];
}



// Returns the densities of ions which have Z index "ir" at radial index ir, for each Z0 (size Z0+1).
const real_t* IonHandler::GetIonDensity(len_t ir, len_t iZ) const{
    real_t *niReturn = new real_t[1+Zs[iZ]];
    const real_t *n_i = unknowns->GetUnknownData(niID);
    for (len_t Z0=0; Z0<Zs[iZ]+1; Z0++)
        niReturn[Z0] = n_i[nr*(ZOffsets[iZ]+Z0) + ir];

    return niReturn;
}

// Returns the total density of ions which have Z index "ir" at radial index ir (summed over Z0).
const real_t IonHandler::GetTotalIonDensity(len_t ir, len_t iZ) const{
    real_t niReturn = 0;
    const real_t *n_i = unknowns->GetUnknownData(niID);
    for (len_t Z0=0; Z0<Zs[iZ]+1; Z0++)
        niReturn += n_i[nr*(ZOffsets[iZ]+Z0) + ir];

    return niReturn;
}


const real_t IonHandler::GetTritiumDensity(len_t ir, len_t *tritiumIndices, len_t numTritiumIndices) const{
    const real_t *n_i = unknowns->GetUnknownData(niID);
    real_t nT;
    if (!(numTritiumIndices==0))
        for (len_t it = 0; it<numTritiumIndices; it++)
            nT += n_i[nr*ZOffsets[tritiumIndices[it]] + ir] + n_i[nr*(1+ZOffsets[tritiumIndices[it]]) + ir];

    return nT; 

}


// Returns the density of ions which are characterised by 
// atomic number Z and charge number Z0 at radial index ir.
const real_t IonHandler::GetFreePlusBoundElectronDensity(len_t ir) const{
    real_t ntot = 0;
    const real_t *n_i = unknowns->GetUnknownData(niID);
    for (len_t iz=0; iz<nZ; iz++)
        for (len_t Z0=0; Z0<Zs[iz]+1; Z0++)
            ntot += Zs[iz]*n_i[nr*(ZOffsets[iz]+Z0) + ir];
    

    return ntot;
}


const real_t IonHandler::GetZeff(len_t ir) const{
    real_t nfreeZ0, nfree;

    const real_t *n_i = unknowns->GetUnknownData(niID);
    for (len_t iz=0; iz<nZ; iz++){
        for (len_t Z0=0; Z0<Zs[iz]+1; Z0++){
            nfree   += Z0*n_i[nr*(ZOffsets[iz]+Z0) + ir];
            nfreeZ0 += Z0*Z0*n_i[nr*(ZOffsets[iz]+Z0) + ir];
        }
    }
    return nfreeZ0/nfree;
}


const real_t IonHandler::GetZtot(len_t ir) const{
    real_t ntotZ, ntot;

    const real_t *n_i = unknowns->GetUnknownData(niID);
    for (len_t iz=0; iz<nZ; iz++){
        for (len_t Z0=0; Z0<Zs[iz]+1; Z0++){
            ntot  += Zs[iz]*n_i[nr*(ZOffsets[iz]+Z0) + ir];
            ntotZ += Zs[iz]*Zs[iz]*n_i[nr*(ZOffsets[iz]+Z0) + ir];
        }
    }
    return ntotZ/ntot;
}





void IonHandler::DeallocateAll(){
    delete [] ZOffsets;
//    delete [] ZList;
//    delete [] Z0List;
}