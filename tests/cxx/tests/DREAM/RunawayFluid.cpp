/**
 * A test for RunawayFluid in DREAM, that verifies the implementation of 
 * derived collision quantities such as the avalanche growth rate and critical E field. 
 * Calculations are benchmarked with values tabulated by DREAM simulations in
 * commit c8f1923d962b3b565ace4e2b033e37ad0a0cb5a8.
 * The Eceff calculation was compared with the function used to generate figures (2-3) 
 * of Hesslow et al, PPCF 60, 074010 (2018), CODE_screened/getEceffWithSynch.m, yielding
 * errors <1% in all three cases when the same bremsstrahlung formula was used in DREAM.
 * The pc calculation was compared with the function behind figure 1 of
 * Hesslow et al Nucl Fusion 59 084004 (2019), which can be found under Linnea's folder
 * CODE_screened/screened_avalanche_implementation/calculateScreenedAvaGrowth
 * and saw agreement within ~1% for all cases tried.
 */

#include <vector>
#include <string>
#include "RunawayFluid.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include <iostream>

using namespace DREAMTESTS::_DREAM;
using namespace std;




/**
 * Run this test.
 */
bool RunawayFluid::Run(bool) {

    bool success = true;
    if (CompareEceffWithTabulated())
        this->PrintOK("The calculation of Eceff agrees with tabulated values.");
    else {
        success = false;
        this->PrintError("The Eceff calculation test failed.");
    }
    if (CompareGammaAvaWithTabulated())
        this->PrintOK("The avalanche growth rate calculation agrees with tabulated values.");
    else {
        success = false;
        this->PrintError("The avalanche growth rate calculation test failed.");
    }

    return success;
}


/**
 * Generate an unknown quantity handler.
 */
DREAM::FVM::UnknownQuantityHandler *RunawayFluid::GetUnknownHandler(DREAM::FVM::Grid *g, 
    const len_t N_IONS, const len_t *Z_IONS, const real_t ION_DENSITY_REF, const real_t T_cold) {
    DREAM::FVM::UnknownQuantityHandler *uqh = new DREAM::FVM::UnknownQuantityHandler();

    len_t nZ0 = 0;
    for (len_t i = 0; i < N_IONS; i++)
        nZ0 += Z_IONS[i] + 1;

    this->id_ions = uqh->InsertUnknown(DREAM::OptionConstants::UQTY_ION_SPECIES, g, nZ0);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_COLD, g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_HOT, g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_TOT, g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_T_COLD, g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_F_HOT, g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_E_FIELD, g);
    

    real_t ni;
    // Set initial values
    const len_t N = nZ0*g->GetNr();
    real_t *nions = new real_t[N];
    real_t ncold = 0;
    real_t ntot = 0;
    len_t ionOffset = 0, rOffset = 0;
    for (len_t iIon = 0; iIon < N_IONS; iIon++) 
        for (len_t Z0 = 0; Z0 <= Z_IONS[iIon]; Z0++, ionOffset++){ 
            ni = (ionOffset+1) * ION_DENSITY_REF;
            ncold += Z0*ni;
            ntot  += Z_IONS[iIon]*ni;
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
    SETVAL(DREAM::OptionConstants::UQTY_N_HOT,  ncold*1e-12);
    SETVAL(DREAM::OptionConstants::UQTY_N_TOT,  ntot);
    SETVAL(DREAM::OptionConstants::UQTY_T_COLD, T_cold);
    SETVAL(DREAM::OptionConstants::UQTY_F_HOT, 0);

    for(len_t i=0;i<g->GetNr();i++)
        temp[i] = 20*(30*i+1);
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_E_FIELD,temp);
//    SETVAL(DREAM::OptionConstants::UQTY_E_FIELD, 0);

    delete [] nions;
    delete [] temp;
    
    return uqh;
}


/**
 * Generate a default ion handler.
 */
DREAM::IonHandler *RunawayFluid::GetIonHandler(
    DREAM::FVM::Grid *g, DREAM::FVM::UnknownQuantityHandler *uqh, const len_t N_IONS, const len_t *Z_IONS
) {
    vector<string> names(N_IONS);
    for (len_t i = 0; i < N_IONS; i++)
        names[i] = "";//ION_NAMES[i];

    return new DREAM::IonHandler(
        g->GetRadialGrid(), uqh, Z_IONS, N_IONS, names
    );
}


DREAM::RunawayFluid *RunawayFluid::GetRunawayFluid(DREAM::CollisionQuantity::collqty_settings *cq,const len_t N_IONS,const len_t *Z_IONS, const real_t ION_DENSITY_REF, const real_t T_cold, const real_t B0, const len_t nr){
    DREAM::FVM::Grid *grid = this->InitializeFluidGrid(nr,B0);
    
    DREAM::FVM::UnknownQuantityHandler *unknowns = GetUnknownHandler(grid,N_IONS, Z_IONS, ION_DENSITY_REF,T_cold);
    DREAM::IonHandler *ionHandler = GetIonHandler(grid,unknowns, N_IONS, Z_IONS);
    DREAM::OptionConstants::momentumgrid_type gridtype = DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI;

    DREAM::CoulombLogarithm *lnLEE = new DREAM::CoulombLogarithm(grid,unknowns,ionHandler,gridtype,cq,DREAM::CollisionQuantity::LNLAMBDATYPE_EE);
    DREAM::CoulombLogarithm *lnLEI = new DREAM::CoulombLogarithm(grid,unknowns,ionHandler,gridtype,cq,DREAM::CollisionQuantity::LNLAMBDATYPE_EI);
    DREAM::SlowingDownFrequency *nuS = new DREAM::SlowingDownFrequency(grid,unknowns,ionHandler,lnLEE,lnLEI,gridtype,cq);
    DREAM::PitchScatterFrequency *nuD = new DREAM::PitchScatterFrequency(grid,unknowns,ionHandler,lnLEI,lnLEE,gridtype,cq);

    DREAM::RunawayFluid *REFluid = new DREAM::RunawayFluid(grid, unknowns, nuS,nuD,lnLEE,lnLEI, cq);    
    REFluid->Rebuild(false);
    return REFluid;
}

bool RunawayFluid::CompareEceffWithTabulated(){

    DREAM::CollisionQuantity::collqty_settings *cq =
        new DREAM::CollisionQuantity::collqty_settings;

    cq->collfreq_type = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED;
    cq->collfreq_mode = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
    cq->lnL_type      = DREAM::OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;
    cq->bremsstrahlung_mode = DREAM::OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_STOPPING_POWER;

    len_t nr = 1;

    real_t Eceff1, Eceff2, Eceff3;
    const len_t N_IONS = 2;
    const len_t Z_IONS[N_IONS] = {10,18};
    real_t ION_DENSITY_REF = 1e18; // m-3
    real_t T_cold = 1; // eV
    real_t B0 = 5;
    DREAM::RunawayFluid *REFluid = GetRunawayFluid(cq,N_IONS, Z_IONS, ION_DENSITY_REF, T_cold,B0,nr);
    Eceff1 = REFluid->GetEffectiveCriticalField(0);

    B0 = 0.1;
    REFluid = GetRunawayFluid(cq,N_IONS, Z_IONS, ION_DENSITY_REF, T_cold,B0,nr);
    Eceff2 = REFluid->GetEffectiveCriticalField(0);

    const len_t N_IONS2 = 1;
    const len_t Z_IONS2[N_IONS2] = {2};
    ION_DENSITY_REF = 1e20; // m-3
    T_cold = 50; // eV
    B0 = 3;
    REFluid = GetRunawayFluid(cq,N_IONS, Z_IONS2, ION_DENSITY_REF, T_cold,B0,nr);
    Eceff3 = REFluid->GetEffectiveCriticalField(0);

    real_t TabulatedEceff1 = 8.88124;
    real_t TabulatedEceff2 = 8.00712;
    real_t TabulatedEceff3 = 91.5226;
    real_t delta1 = abs(Eceff1-TabulatedEceff1)/TabulatedEceff1;
    real_t delta2 = abs(Eceff2-TabulatedEceff2)/TabulatedEceff2;
    real_t delta3 = abs(Eceff3-TabulatedEceff3)/TabulatedEceff3;

    real_t threshold = 1e-2;

/*
    cout << "Delta1: " << delta1 << endl;
    cout << "Delta2: " << delta2 << endl;
    cout << "Delta3: " << delta3 << endl;
*/  
    return (delta1 < threshold) && (delta2 < threshold) && (delta3 < threshold);
}


/**
 * Evalutes the semi-analytic avalanche growth rate in a Neon-Argon plasma for 
 * three different E fields and compares with tabulated values.
 */
bool RunawayFluid::CompareGammaAvaWithTabulated(){

    DREAM::CollisionQuantity::collqty_settings *cq =
        new DREAM::CollisionQuantity::collqty_settings;

    cq->collfreq_type = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED;
    cq->collfreq_mode = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
    cq->lnL_type      = DREAM::OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;
    cq->bremsstrahlung_mode = DREAM::OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT;
    cq->pstar_mode = DREAM::OptionConstants::COLLQTY_PSTAR_MODE_COLLISIONLESS;

    len_t nr = 3;
    const len_t N_IONS = 2;
    const len_t Z_IONS[N_IONS] = {10,18};
    real_t ION_DENSITY_REF = 1e18; // m-3
    real_t T_cold = 1; // eV
    real_t B0 = 5;
    DREAM::RunawayFluid *REFluid = GetRunawayFluid(cq,N_IONS, Z_IONS, ION_DENSITY_REF, T_cold,B0,nr);

    const real_t *GammaAva =  REFluid->GetAvalancheGrowthRate();
    const real_t GammaTabulated[3] = {161.106, 11778.7, 25054.8};
    real_t *deltas = new real_t[3];
    for(len_t ir=0; ir<nr;ir++)
        deltas[ir] = abs(GammaAva[ir]-GammaTabulated[ir])/GammaTabulated[ir];

/*
    cout << "Delta1: " << deltas[1] << endl;
    cout << "Delta2: " << deltas[2] << endl;
    cout << "Delta3: " << deltas[3] << endl;
*/
    real_t threshold = 2e-2;
    return (deltas[1]<threshold) && (deltas[2]<threshold) && (deltas[3]<threshold);    
}