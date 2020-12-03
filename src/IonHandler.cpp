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
#include "DREAM/Constants.hpp"

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


// Standard atomic weight of the elements, which are the atomic masses weighted
// by natural abundances and normalized to the atomic mass unit (dalton) u
const real_t IonHandler::atomicMassInMu[nIonMass] = 
{
// Z = 1,      2,      3,      4,      5,      6,      7,      8,      9,      10 
       1.008,  4.0026, 6.94,   9.0122, 10.81,  12.011, 14.007, 15.999, 18.998, 20.180, 
//     11,     12,     13,     14,     15,     16,     17,     18,     19,     20,        
       22.990, 24.305, 26.982, 28.085, 30.974, 32.06,  35.45,  39.95,  39.098, 40.078, 
//     21,     22,     23,     24,     25,     26,     27,     28,     29,     30   
       44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 
//     31,     32,     33,     34,     35,     36,     37,     38,     39,     40
       69.723, 72.63,  74.922, 78.971, 79.904, 83.798, 85.468, 87.62,  88.906, 91.224
};

/**
 * Constructor.
 *
 * rg:      Radial grid on which the ions live.
 * u:       List of unknown quantities.
 * Z:       List of atomic charges for each species (size NZ).
 * NZ:      Number of atomic species.
 * names:   List of strings defining the names by which each ion species
 *          will be referred (must have 'NZ' elements).
 * tritium: List of names of the tritium ion species.
 */
IonHandler::IonHandler(
    FVM::RadialGrid *rg, FVM::UnknownQuantityHandler *u, const len_t *Z, len_t NZ,
    vector<string>& names, vector<string>& tritium
) {
    rGrid = rg;
    unknowns = u;

    this->Zs  = Z;
    nr = rGrid->GetNr();
    nZ = NZ;

    niID = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

    this->ionNames = names;
    this->tritiumNames = tritium;
    this->nTritium = tritium.size();
    this->tritiumIndices = new len_t[nTritium];

    // Find index of tritum ions
    len_t ti = 0;
    for (len_t t = 0; t < tritium.size(); t++){
        for (len_t i = 0; i < names.size(); i++) 
            if (tritium[t] == names[i]) {
                this->tritiumIndices[ti++] = i;
                break;
            }
        if (ti != t+1)
            throw FVM::FVMException("Species '%s' declared as tritium, but ion species has not been defined.", tritium[t].c_str());
    }
    Initialize();
}

/**
 * Destructor.
 */
IonHandler::~IonHandler(){
    DeallocateAll();
    delete [] tritiumIndices;
}


/**
 * Calculate number of ion charge states and determine the indices
 * for where different ion species start in the ion array.
 */
void IonHandler::Initialize() {
    DeallocateAll();
    nzs = 0;
    for (len_t it=0; it<nZ; it++)
        nzs += Zs[it]+1;

    ZOffsets = new len_t[nZ];
    ZOffsets[0] = 0;
    if (nZ>0)
        for (len_t iz=1; iz<nZ; iz++)
            ZOffsets[iz] = ZOffsets[iz-1] + Zs[iz-1] + 1;

    izsList = new len_t[nzs];
    Z0sList = new len_t[nzs];
    for(len_t iz=0; iz<nZ; iz++)
        for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
            len_t indZ = GetIndex(iz,Z0);
            izsList[indZ] = iz;
            Z0sList[indZ] = Z0;
        }

    nfree  = new real_t[nr];
    ntot   = new real_t[nr];
    nbound = new real_t[nr];
    nZ0Z0  = new real_t[nr];
    nZZ    = new real_t[nr];
    Zeff   = new real_t[nr];
    Ztot   = new real_t[nr];

    mi = new real_t[nZ];
    for(len_t iz=0; iz<nZ; iz++){
        if(Zs[iz]==1){ // assume pure deuterium unless it is marked as tritium
            bool isTritium = false;
            for(len_t it=0; it<nTritium; it++)
                if(iz==tritiumIndices[it])
                    isTritium = true;
            mi[iz] = isTritium ? Constants::mT : Constants::mD;
        } else if ( Zs[iz] > nIonMass ) // if heavier species than we store data for, assume simple linear scaling
            mi[iz] = 2.3*Zs[iz] * Constants::mu; 
        else // read from table
            mi[iz] = atomicMassInMu[Zs[iz]-1] * Constants::mu;
    }
}


/**
 * Evaluates and stores various moments of the ion densities. 
 * For initialization, this should be called immediately after
 * ion densities have been initialized. 
 */
void IonHandler::Rebuild(){
    for(len_t ir=0; ir<nr; ir++){
        nfree[ir]  = 0;
        ntot[ir]   = 0;
        nbound[ir] = 0;
        Zeff[ir]   = 0;
        Ztot[ir]   = 0;
        nZ0Z0[ir]  = 0;
        nZZ[ir]    = 0;

        for (len_t iz = 0; iz < nZ; iz++)
            for (len_t Z0 = 0; Z0<=Zs[iz]; Z0++){
                real_t ni = GetIonDensity(ir,iz,Z0);
                nfree[ir]  += Z0*ni;
                nbound[ir] += (Zs[iz]-Z0)*ni;
                ntot[ir]   += Zs[iz]*ni;
                nZ0Z0[ir]  += Z0*Z0*ni;
                nZZ[ir]    += Zs[iz]*Zs[iz]*ni;
            }
        if(nfree[ir])
            Zeff[ir] = nZ0Z0[ir] / nfree[ir];
        if(ntot[ir])
            Ztot[ir] = nZZ[ir] / ntot[ir];    

    }
}


/**
 *  Returns the density of ions which are characterised by 
 * atomic number Z and charge number Z0 at radial index ir.
 */
const real_t IonHandler::GetIonDensityAtZ(len_t ir, len_t Z, len_t Z0) const{
    real_t niReturn = 0;
    const real_t *n_i = unknowns->GetUnknownData(niID);
    for (len_t iz=0; iz<nZ; iz++)
        if (Zs[iz] == Z){
            len_t Zind = GetIndex(iz,Z0);
            niReturn += n_i[nr*Zind + ir];
        }
    
    return niReturn;
}

// Returns the density of ions which are characterised by 
// atomic number Z and charge number Z0 at radial index ir.
const real_t IonHandler::GetIonDensity(len_t ir, len_t iz, len_t Z0) const{

    if (Z0 > Zs[iz])
        throw FVM::FVMException("Ion charge number cannot be larger than atomic number.");

    const real_t *n_i = unknowns->GetUnknownData(niID);
    len_t Zind = GetIndex(iz,Z0);
    return n_i[nr*Zind + ir];
}



// Returns the densities of ions which have Z index "ir" at radial index ir, for each Z0 (size Z0+1).
const real_t* IonHandler::GetIonDensity(len_t ir, len_t iZ) const{
    real_t *niReturn = new real_t[1+Zs[iZ]];
    const real_t *n_i = unknowns->GetUnknownData(niID);
    for (len_t Z0=0; Z0<Zs[iZ]+1; Z0++){
        len_t Zind = GetIndex(iZ,Z0);
        niReturn[Z0] = n_i[nr*Zind + ir];
    }

    return niReturn;
}

// Returns the total density of ions which have Z index "ir" at radial index ir (summed over Z0).
const real_t IonHandler::GetTotalIonDensity(len_t ir, len_t iZ) const{
    real_t niReturn = 0;
    const real_t *n_i = unknowns->GetUnknownData(niID);
    for (len_t Z0=0; Z0<Zs[iZ]+1; Z0++){
        len_t Zind = GetIndex(iZ,Z0);
        niReturn += n_i[nr*Zind + ir];
    }
    return niReturn;
}


/**
 * Calculates the density of tritium in the plasma.
 *
 * ir: Radius at which to calculate the tritium density.
 */
const real_t IonHandler::GetTritiumDensity(len_t ir) const {
    const real_t *n_i = unknowns->GetUnknownData(niID);

    real_t nT = 0;
    if (this->nTritium > 0) {
        for (len_t it = 0; it < this->nTritium; it++)
            nT += n_i[nr*ZOffsets[this->tritiumIndices[it]] + ir] +     // Z0 = 0
                  n_i[nr*(1+ZOffsets[this->tritiumIndices[it]]) + ir];  // Z0 = 1
    }

    return nT; 
}


/**
 * The inverse of GetIndex(...): takes the ion index and 
 * returns the corresponding iz and Z0
 */
void IonHandler::GetIonIndices(len_t nMultiple, len_t &iz_in, len_t &Z0_in){
    iz_in = izsList[nMultiple];
    Z0_in = Z0sList[nMultiple];
}


/**
 * Checks whether the ion with the given index is a tritium
 * ion species.
 */
bool IonHandler::IsTritium(const len_t iIon) const {
    for (len_t i = 0; i < this->nTritium; i++) {
        if (iIon == this->tritiumIndices[i])
            return true;
    }

    return false;
}


/**
 * Calculate the effective bound charge,
 *   Zeff0 = sum_i( n_i*(Z_i^2 - Z_i0^2) ) / ntot
 * where
 *   ntot = sum_i( n_i*Z_i )
 * and the indices run over all ion species and charge states.
 */
real_t *IonHandler::evaluateZeff0() {
    real_t *Zeff0 = new real_t[nr];

    for (len_t ir = 0; ir < nr; ir++)
        Zeff0[ir] = evaluateZeff0(ir);

    return Zeff0;
}


/**
 * Calculate the effective bound charge.
 */
real_t IonHandler::evaluateZeff0(len_t ir) {
    real_t ntot   = 0;
    real_t ntotZ0 = 0;
    for (len_t iz = 0; iz < nZ; iz++) {
        real_t ntotz = GetTotalIonDensity(ir, iz);
        ntot += Zs[iz]*ntotz;

        for (len_t Z0=1; Z0<Zs[iz]+1; Z0++)
            ntotZ0 += (Zs[iz]*Zs[iz] - Z0*Z0)*ntotz;
    }
    return ntotZ0 / ntot;
}


/**
 * Evaluate
 *   Z0Z = sum_i( n_i*Z_i0*Z_i ) / ntot
 * where
 *   ntot = sum_i( n_i*Z_i )
 * and the indices run over all ion species and charge states.
 */
real_t *IonHandler::evaluateZ0Z() {
    real_t *Z0Z = new real_t[nr];
    for (len_t ir = 0; ir < nr; ir++)
        Z0Z[ir] = evaluateZ0Z(ir);
    return Z0Z;
}


/**
 * Evaluate Z0Z at the given radius.
 */
real_t IonHandler::evaluateZ0Z(len_t ir) {
    real_t ntot   = 0;
    real_t nZmul  = 0;
    for (len_t iz = 0; iz < nZ; iz++) {
        real_t ntotz = GetTotalIonDensity(ir, iz);
        ntot += Zs[iz]*ntotz;
        for (len_t Z0=1; Z0<Zs[iz]+1; Z0++)
            nZmul += (Z0 * Zs[iz])*ntotz;
    }
    return nZmul / ntot;
}

/**
 * Evaluate
 *   Z0_Z = sum_i( n_i*Z_i0 / Z_i ) / ntot
 * where
 *   ntot = sum_i( n_i*Z_i )
 * and the indices run over all ion species and charge states.
 */
real_t *IonHandler::evaluateZ0_Z() {
    real_t *Z0_Z = new real_t[nr];
    for (len_t ir = 0; ir < nr; ir++)
        Z0_Z[ir] = evaluateZ0_Z(ir);
    return Z0_Z;
}


/**
 * Evaluate Z0_Z at the specified radius.
 */
real_t IonHandler::evaluateZ0_Z(len_t ir) {
    real_t ntot   = 0;
    real_t nZdiv  = 0;
    for (len_t iz = 0; iz < nZ; iz++) {
        real_t ntotz = GetTotalIonDensity(ir, iz);
        ntot += Zs[iz]*ntotz;
        for (len_t Z0=1; Z0<Zs[iz]+1; Z0++)
            nZdiv += (Z0 / Zs[iz])*ntotz;
    }
    return nZdiv / ntot;
}


void IonHandler::DeallocateAll(){
    if(ZOffsets == nullptr)
        return;
    delete [] ZOffsets;
    delete [] izsList;
    delete [] Z0sList;
    delete [] nfree;
    delete [] ntot;
    delete [] nbound;
    delete [] nZ0Z0;
    delete [] nZZ;
    delete [] Zeff;
    delete [] Ztot;
    delete [] mi;
}
