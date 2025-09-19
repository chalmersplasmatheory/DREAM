/**
 * A test for the BootstrapCurrent term in DREAM.
 */

#include <vector>
#include <string>
#include "BootstrapCurrent.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/BootstrapCurrent.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CoulombLogarithm.hpp"


using namespace DREAMTESTS::_DREAM;
using namespace std;

// TODO: Look at the IonRateEquation for inspiration, 
const len_t N_IONS = /*TODO*/; // Number of ion species
const len_t Z_IONS[N_IONS] = {/*TODO*/;}; // Vector with atomic number for each ion species, # elements == N_IONS
const char ION_NAMES[N_IONS][3] = {/*TODO*/}; // Vector with name of each ion species, # elements == N_IONS


/**
 * Generate a default ion handler.
 */
DREAM::IonHandler *BootstrapCurrent::GetIonHandler(
    DREAM::FVM::Grid *g, DREAM::FVM::UnknownQuantityHandler *uqh
) {
    vector<string> tritiumNames(0), hydrogenNames(0); // TODO: not sure what to do here, maybe need to change this
    vector<string> names(N_IONS);
    len_t *Z = new len_t[N_IONS];// The ion charge numbers must be provided to the IONHandler as a dynamically allocated array to avoid memory issues
    for (len_t i = 0; i < N_IONS; i++){
        names[i] = ION_NAMES[i];
        Z[i]=Z_IONS[i];
    }

    return new DREAM::IonHandler(
        g->GetRadialGrid(), uqh, Z, N_IONS, names, tritiumNames, hydrogenNames
    );
}

/**
 * Generate an unknown quantity handler.
 */
DREAM::FVM::UnknownQuantityHandler *BootstrapCurrent::GetUnknownHandler(DREAM::FVM::Grid *g, bool withIonEnergy) {
    DREAM::FVM::UnknownQuantityHandler *uqh = new DREAM::FVM::UnknownQuantityHandler();

    len_t nZ0 = 0;
    for (len_t i = 0; i < N_IONS; i++)
        nZ0 += Z_IONS[i] + 1;

    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_ION_SPECIES, "0", g, nZ0);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_COLD, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_T_COLD, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_NI_DENS, "0", g, N_IONS);
    if (withIonEnergy)
        uqh->InsertUnknown(DREAM::OptionConstants::UQTY_WI_ENER, "0", g, N_IONS);


    // Set initial values
    const len_t Nsum = N_IONS*g->GetNr();
    const len_t N = nZ0*g->GetNr();;
    real_t *Nions = new real_t[Nsum];
    real_t *nions = new real_t[N];
    len_t  rOffset = 0;
    for (len_t iIon = 0; iIon < N_IONS; iIon++) {
        for (len_t Z0 = 0; Z0 <= Z_IONS[iIon]; Z0++) {
            for (len_t ir = 0; ir < g->GetNr(); ir++, rOffset++){
                nions[rOffset] = /*TODO*/;
                Nions[iIon*nr+ir] += nions[rOffset];
            }
                
        }
    }
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_ION_SPECIES, nions);
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_NI_DENS, nions);

    if (withIonEnergy){
        // TODO: Set ion energies
        real_t *wions = new real_t[Nsum];
        len_t rOffset = 0;
        for (len_t iIon = 0; iIon < N_IONS; iIon++) {
            for (len_t ir = 0; ir < g->GetNr(); ir++, rOffset++)
                wions[rOffset] = /*TODO*/;
        }
        uqh->SetInitialValue(DREAM::OptionConstants::UQTY_WI_ENER, wions);
    }

    // TODO: Set electron densities and energies
    // Set electron quantities
    real_t *ncold = new real_t[g->GetNr()];
    real_t *Tcold = new real_t[g->GetNr()];
    for (len_t i = 0; i < g->GetNr(); i++){
        ncold[i] = /*TODO*/;
        Tcold[i] = /*TODO*/;
    }
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_N_COLD, /*TODO*/);
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_T_COLD, /*TODO*/);

    delete [] nions;
    delete [] Nions;
    if (withIonEnergy)
        delete [] wions;
    delete [] ncold;
    delete [] Tcold;
    
    return uqh;
}

/**
 * Check if the bootstrap current is consistent with IDA.
 */
bool BootstrapCurrent::CheckBootstrap(bool withIonEnergy) {
    real_t successRelErrorThreshold = 1e-5; // TODO: Probably want to change this
    bool success = true;
    const len_t Nr = /*TODO*/; 
    DREAM::FVM::Grid *grid = this->InitializeFluidGrid(Nr);
    DREAM::FVM::UnknownQuantityHandler *uqh = GetUnknownHandler(grid, withIonEnergy);
    DREAM::IonHandler *ih = GetIonHandler(grid, uqh);
    ih->Rebuild();

    DREAM::OptionConstants::momentumgrid_type gridtype = DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI;
    DREAM::CollisionQuantity::collqty_settings *cq = new DREAM::CollisionQuantity::collqty_settings;

    cq->collfreq_type = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED;
    cq->collfreq_mode = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
    cq->lnL_type      = DREAM::OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;
    cq->bremsstrahlung_mode = DREAM::OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_STOPPING_POWER;

    DREAM::CoulombLogarithm *lnLEE = new DREAM::CoulombLogarithm(grid,uqh,ih,gridtype,cq,DREAM::CollisionQuantity::LNLAMBDATYPE_EE);
    DREAM::BootstrapCurrent *bootstrap = new DREAM::BootstrapCurrent(grid, uqh, ih, lnLEE);
    bootstrap->Rebuild();

    // Construct solution vector
    real_t *j_bs = new real_t[Nr];
    real_t *j_bs_IDA = new real_t[Nr];
    
    // Clear 'j_bs'
    for (len_t ir = 0; ir < Nr; ir++){
        j_bs[ir] = 0;
        j_bs_IDA[ir] = /*TODO*/;
    }
    
    // Construct equation for each ion species
    DREAM::BootstrapElectronDensityTerm term_ncold = new DREAM::BootstrapElectronDensityTerm(grid, uqh, bootstrap, ih);
    DREAM::BootstrapElectronTemperatureTerm term_Tcold = new DREAM::BootstrapElectronTemperatureTerm(grid, uqh, bootstrap, ih);

    DREAM::BootstrapIonDensityTerm *term_Ni[N_IONS];
    if (withIonEnergy)
        DREAM::BootstrapIonThermalEnergyTerm *term_wi[N_IONS];
    for (len_t iIon = 0; iIon < N_IONS; iIon++){
        term_Ni[iIon] = new DREAM::BootstrapIonDensityTerm(grid, uqh, bootstrap, ih, iIon);
        if (withIonEnergy)
            term_wi[iIon] = new DREAM::BootstrapIonThermalEnergyTerm(grid, uqh, bootstrap, ih, iIon);
    }

    
    term_ncold->Rebuild(nullptr, nullptr, uqh);
    term_ncold->SetVectorElements(j_bs, nullptr);
    term_Tcold->Rebuild(nullptr, nullptr, uqh);
    term_Tcold->SetVectorElements(j_bs, nullptr);
    len_t rOffset = 0;
    for (len_t iIon = 0; iIon < N_IONS; iIon++) {
        term_Ni[iIon]->Rebuild(nullptr, nullptr, uqh);
        term_Ni[iIon]->SetVectorElements(j_bs, nullptr);
        if (withIonEnergy){
            term_wi[iIon]->Rebuild(nullptr, nullptr, uqh);
            term_wi[iIon]->SetVectorElements(j_bs, nullptr);
        }
    }

    /* 
     * TODO: Alternatively, the IDA data for bootstrap current can be put in the vector j_bs_IDA in a loop here
     */

    // Sum of all elements should vanish
    real_t deltas;
    for (len_t ir = 0; ir < Nr; ir++) {
        deltas = abs(j_bs[ir] - j_bs_IDA[ir]);
        if(deltas>successRelErrorThreshold){
            success = false;
            cout << "delta [rel error]: " << deltas << endl;
        }
    }


    // Deallocate equation terms
    delete term_ncold;
    delete term_Tcold;
    for (len_t iIon = 0; iIon < N_IONS; iIon++){
        delete term_Ni[iIon];
        if (withIonEnergy)
            delete term_wi[iIon];
    }
    delete [] j_bs;
    delete [] j_bs_IDA;
    delete bootstrap;
    delete lnLEE;
    delete cq;
    delete ih;
    delete uqh;
    delete grid;

    return success;
}

/**
 * Run this test.
 */
bool BootstrapCurrent::Run(bool) {
    bool success = true;

    if ((success &= CheckBootstrap(true)) && (success &= CheckBootstrap(false)))
        this->PrintOK("The bootstrap current is consistent with IDA.");
    else
        this->PrintError("The bootstrap current is not consistent with IDA.");

    return success;
}

