/**
 * Implementation of the object which helps map ion densities to
 * ion names and charges. This is needed since the ion densities
 * are stored a single, monolithic density in the EquationSystem.
 * Since the UnknownQuantity doesn't contain any information about
 * the ion charge, we must keep that information in this separate
 * object instead.
 */

#include <vector>
#include <string>
#include "DREAM/IonHandler.hpp"

using namespace DREAM;
using namespace std;


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


/**
 * Constructor.
 *
 * rg:    Radial grid on which the ions live.
 * u:     List of unknown quantities.
 * Z:     List of atomic charges for each species (size NZ).
 * NZ:    Number of atomic species.
 * names: List strings defining the names by which each ion species
 *        will be referred (must have 'NZ' elements).
 */
IonHandler::IonHandler(
    FVM::RadialGrid *rg, FVM::UnknownQuantityHandler *u, const len_t *Z, len_t NZ,
    vector<string>& names
) {
    rGrid = rg;
    unknowns = u;

    this->Zs  = Z;
    nr = rGrid->GetNr();
    nZ = NZ;

    niID = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

    this->ionNames = names;

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

    if (Z0 > Zs[iz]){
        throw FVM::FVMException("Ion charge number cannot be larger than atomic number.");
    }

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
    if (numTritiumIndices>0)
        for (len_t it = 0; it<numTritiumIndices; it++)
            nT += n_i[nr*ZOffsets[tritiumIndices[it]] + ir] + n_i[nr*(1+ZOffsets[tritiumIndices[it]]) + ir];

    return nT; 

}


/**
 * Returns the total electron density n_tot = n_free + n_bound.
 *
 * ntot: If NOT 'nullptr', this array contains the total electron
 *       density upon return. Otherwise, if 'nullptr', new memory
 *       is allocated for n_tot.
 */
real_t* IonHandler::evaluateFreePlusBoundElectronDensityFromQuasiNeutrality(real_t *ntot){
    if (ntot == nullptr)
        ntot = new real_t[nr];

    // Initialize array to zero
    for (len_t ir = 0; ir < nr; ir++)
        ntot[ir] = 0;

    for (len_t ir = 0; ir < nr; ir++)
        for (len_t iz = 0; iz < nZ; iz++)
            ntot[ir] += Zs[iz]*GetTotalIonDensity(ir,iz);

    return ntot;
}



/**
 * Returns the free electron density n_free.
 *
 * nfree: If NOT 'nullptr', this array contains the free electron
 *       density upon return. Otherwise, if 'nullptr', new memory
 *       is allocated for n_free.
 */
real_t* IonHandler::evaluateFreeElectronDensityFromQuasiNeutrality(real_t *nfree){
    if (nfree == nullptr)
        nfree = new real_t[nr];

    // Initialize array to zero
    for (len_t ir = 0; ir < nr; ir++)
        nfree[ir] = 0;

    for (len_t ir = 0; ir < nr; ir++)
        for (len_t iz = 0; iz < nZ; iz++)
            for (len_t Z0 = 1; Z0<=Zs[iz]; Z0++)
                nfree[ir] += Z0*GetIonDensity(ir,iz,Z0);

    return nfree;
}



/**
 * Returns the bound electron density n_bound.
 *
 * nbound: If NOT 'nullptr', this array contains the free electron
 *         density upon return. Otherwise, if 'nullptr', new memory
 *         is allocated for n_free.
 */
real_t* IonHandler::evaluateBoundElectronDensityFromQuasiNeutrality(real_t *nbound){
    if (nbound == nullptr)
        nbound = new real_t[nr];

    // Initialize array to zero
    for (len_t ir = 0; ir < nr; ir++)
        nbound[ir] = 0;

    for (len_t ir = 0; ir < nr; ir++)
        for (len_t iz = 0; iz < nZ; iz++)
            for (len_t Z0 = 1; Z0<=Zs[iz]; Z0++)
                nbound[ir] += (GetZ(iz) - Z0)*GetIonDensity(ir,iz,Z0);

    return nbound;
}






real_t* IonHandler::evaluateZeff(){
    real_t nfreeZ0, nfree;

    real_t *Zeff = new real_t[nr];
    for(len_t ir=0; ir<nr; ir++){
        for (len_t iz=0; iz<nZ; iz++){
            for (len_t Z0=0; Z0<Zs[iz]+1; Z0++){
                nfree   += Z0*GetIonDensity(ir,iz,Z0);
                nfreeZ0 += Z0*Z0*GetIonDensity(ir,iz,Z0);
            }
        }
        Zeff[ir] = nfreeZ0/nfree;
    }
    return Zeff;
}


real_t* IonHandler::evaluateZtot(){
    real_t ntotZ;
    real_t *ntot = evaluateFreePlusBoundElectronDensityFromQuasiNeutrality();
    real_t *Ztot = new real_t[nr];
    for (len_t ir=0; ir<nr; ir++){
        ntotZ = 0;
        for (len_t iz=0; iz<nZ; iz++){
            ntotZ += Zs[iz]*Zs[iz]*GetTotalIonDensity(ir,iz);
        }
        Ztot[ir] = ntotZ/ntot[ir];
    }
    delete [] ntot;
    return Ztot;
}





void IonHandler::DeallocateAll(){
    delete [] ZOffsets;
}
