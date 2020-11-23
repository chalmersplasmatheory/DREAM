/**
 * A test to verify the updated values for the mean excitation energy,
 * using the model by Sauer, Sabin, Oddershede J Chem Phys 148, 174307 (2018)
 */

#include <vector>
#include <string>
#include "MeanExcitationEnergy.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CoulombLogarithm.hpp"
#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include <iostream>

using namespace DREAMTESTS::_DREAM;
using namespace std;


/**
 * Run this test.
 */
bool MeanExcitationEnergy::Run(bool) {

    bool success = true;
    if (CompareMeanExcitationEnergyWithTabulated())
        this->PrintOK("The mean excitation energy agrees with tabulated values.");
    else {
        success = false;
        this->PrintError("The mean excitation energy test failed.");
    }
    return success;
}


/**
 * Generate an unknown quantity handler.
 */
DREAM::FVM::UnknownQuantityHandler *MeanExcitationEnergy::GetUnknownHandler(DREAM::FVM::Grid *g, 
    const len_t N_IONS, const len_t *Z_IONS, const real_t ION_DENSITY_REF, const real_t T_cold) {
    DREAM::FVM::UnknownQuantityHandler *uqh = new DREAM::FVM::UnknownQuantityHandler();

    len_t nZ0 = 0;
    for (len_t i = 0; i < N_IONS; i++)
        nZ0 += Z_IONS[i] + 1;

    this->id_ions = uqh->InsertUnknown(DREAM::OptionConstants::UQTY_ION_SPECIES, "0", g, nZ0);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_COLD, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_T_COLD, "0", g);
    
    real_t ni;
    // Set initial values
    // we don't care about the density, just set it to the same value everywhere.
    const len_t N = nZ0*g->GetNr();
    real_t *nions = new real_t[N];
    real_t ncold = 0;
    len_t ionOffset = 0, rOffset = 0;
    for (len_t iIon = 0; iIon < N_IONS; iIon++) 
        for (len_t Z0 = 0; Z0 <= Z_IONS[iIon]; Z0++, ionOffset++){ 
            ni =  ION_DENSITY_REF;
            ncold += Z0*ni;
            for (len_t ir = 0; ir < g->GetNr(); ir++, rOffset++)
                nions[rOffset] = ni;
        }
        
    
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_ION_SPECIES, nions);

    #define SETVAL(NAME, v) do { \
            for (len_t i = 0; i < g->GetNr(); i++) \
                temp[i] = (v); \
            uqh->SetInitialValue((NAME), temp); \
        } while (false)

    // Set electron quantities
    real_t *temp = new real_t[g->GetNr()];
    SETVAL(DREAM::OptionConstants::UQTY_N_COLD, ncold);
    SETVAL(DREAM::OptionConstants::UQTY_T_COLD, T_cold);

    delete [] nions;
    delete [] temp;
    
    return uqh;
}


/**
 * Generate a default ion handler.
 */
DREAM::IonHandler *MeanExcitationEnergy::GetIonHandler(
    DREAM::FVM::Grid *g, DREAM::FVM::UnknownQuantityHandler *uqh, const len_t N_IONS, const len_t *Z_IONS
) {
    vector<string> tritiumNames(0);
    vector<string> names(N_IONS);
    for (len_t i = 0; i < N_IONS; i++)
        names[i] = "";

    return new DREAM::IonHandler(
        g->GetRadialGrid(), uqh, Z_IONS, N_IONS, names, tritiumNames
    );
}

void MeanExcitationEnergy::GetMeanExcitationEnergies(real_t *meanExcitationEnergies, DREAM::CollisionQuantity::collqty_settings *cq, const len_t N_IONS, const len_t *Z_IONS, 
    const len_t N_SPECIES_TO_TEST, const len_t *Z_TO_TEST, const len_t *Z0_TO_TEST,
    const real_t ION_DENSITY_REF, const real_t T_cold, const real_t B0, const len_t nr
){
    DREAM::FVM::Grid *grid = this->InitializeFluidGrid(nr,B0);
    
    DREAM::FVM::UnknownQuantityHandler *unknowns = GetUnknownHandler(grid,N_IONS, Z_IONS, ION_DENSITY_REF,T_cold);
    DREAM::IonHandler *ionHandler = GetIonHandler(grid,unknowns, N_IONS, Z_IONS);
    ionHandler->Rebuild();
    DREAM::OptionConstants::momentumgrid_type gridtype = DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI;

    DREAM::CoulombLogarithm lnLEE(grid,unknowns,ionHandler,gridtype,cq,DREAM::CollisionQuantity::LNLAMBDATYPE_EE);
    DREAM::CoulombLogarithm lnLEI(grid,unknowns,ionHandler,gridtype,cq,DREAM::CollisionQuantity::LNLAMBDATYPE_EI);
    DREAM::SlowingDownFrequency nuS(grid,unknowns,ionHandler,&lnLEE,&lnLEI,gridtype,cq);
    nuS.RebuildRadialTerms();
    
    len_t iz = 0;
    for (len_t is = 0; is < N_SPECIES_TO_TEST; is++) {
        if (Z_TO_TEST[is] != Z_IONS[iz]){ iz++; }
        meanExcitationEnergies[is] = nuS.GetMeanExcitationEnergy(iz,Z0_TO_TEST[is]);
    }

    delete ionHandler;
    delete unknowns;
    delete grid;
}

bool MeanExcitationEnergy::CompareMeanExcitationEnergyWithTabulated(){

    DREAM::CollisionQuantity::collqty_settings *cq =
        new DREAM::CollisionQuantity::collqty_settings;

    cq->collfreq_type = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED;
    cq->collfreq_mode = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
    cq->lnL_type      = DREAM::OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;
    cq->bremsstrahlung_mode = DREAM::OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_STOPPING_POWER;

    len_t nr = 1;

    const len_t N_IONS = 3;
    const len_t Z_IONS[N_IONS] = {10,18,36};
    constexpr int_t N_SPECIES_TO_TEST = 10;
    const len_t Z_TO_TEST[N_SPECIES_TO_TEST] =  {10,10,10,18,18,18,36,36,36,36};
    const len_t Z0_TO_TEST[N_SPECIES_TO_TEST] = { 0, 1, 5, 0, 1, 9, 0, 1,18,30};

    const real_t TABULATED_MEAN_EXCITATION_ENERGIES[] = {
        //Ne0,       Ne1,         Ne5,         Ar0,         Ar1,         Ar9,         Kr0,         Kr1,         Kr18,        Kr30
        2.68494e-04, 3.23288e-04, 6.90021e-04, 3.69277e-04, 4.29551e-04, 1.55969e-03, 7.04502e-04, 8.06097e-04, 3.45034e-03, 1.01486e-02
    };

    real_t ION_DENSITY_REF = 1e18; // m-3
    real_t T_cold = 1; // eV
    real_t B0 = 5;

    real_t meanExcitationEnergies[N_SPECIES_TO_TEST];
    GetMeanExcitationEnergies(meanExcitationEnergies, cq,N_IONS, Z_IONS,
    N_SPECIES_TO_TEST, Z_TO_TEST, Z0_TO_TEST,ION_DENSITY_REF, T_cold,B0,nr);
    
    real_t delta;
    bool success = true;
    const real_t TOLERANCE = 1e-4;
    for (len_t is = 0; is < N_SPECIES_TO_TEST; is++) {
        delta = abs(meanExcitationEnergies[is] - TABULATED_MEAN_EXCITATION_ENERGIES[is]) / TABULATED_MEAN_EXCITATION_ENERGIES[is];
        success = success & (delta < TOLERANCE);
    }

    delete cq;
    return success;
}