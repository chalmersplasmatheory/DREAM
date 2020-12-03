/**
 * A test for RunawayFluid in DREAM, that verifies the implementation of 
 * derived collision quantities such as the avalanche growth rate and critical E field. 
 * Calculations are benchmarked with values tabulated by DREAM simulations in
 * commit c8f1923d962b3b565ace4e2b033e37ad0a0cb5a8.
 * The Eceff values were updated in commit b5c5dab98742f3925c71177e27b536a4693a25aa when
 * the mean excitation energies were updated, and in 32228b5214d54ff69033e2dbd848b7d5daa18a01
 * when the iterative solution was implemented. 
 * *
 * The Eceff calculation was compared with the function used to generate figures (2-3) 
 * of Hesslow et al, PPCF 60, 074010 (2018), CODE_screened/getEceffWithSynch.m, yielding
 * errors <1% in all three cases when the same bremsstrahlung formula was used in DREAM.
 * The iterative calculation of Eceff was compared with the MATLAB implementation of Eceff, 
 * as well as figures 2-3 in the paper above, giving the same values to at least 5 decimals.
 * 
 * The pc calculation was compared with the function behind figure 1 of
 * Hesslow et al Nucl Fusion 59 084004 (2019), which can be found under Linnea's folder
 * CODE_screened/screened_avalanche_implementation/calculateScreenedAvaGrowth
 * and saw agreement within ~1% for all cases tried.
 */

#include <vector>
#include <string>
#include "RunawayFluid.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/ConnorHastie.hpp"
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

    if (CompareConnorHastieRateWithTabulated())
        this->PrintOK("The Connor-Hastie runaway rate is calculated correctly.");
    else {
        success = false;
        this->PrintError("The Connor-Hastie runaway rate test failed.");
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

    this->id_ions = uqh->InsertUnknown(DREAM::OptionConstants::UQTY_ION_SPECIES, "0", g, nZ0);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_COLD, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_HOT, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_TOT, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_T_COLD, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_F_HOT, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_E_FIELD, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_J_TOT, "0", g);
    

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

    delete [] nions;
    delete [] temp;
    
    return uqh;
}

/**
 * Generate an unknown quantity handler.
 */
DREAM::FVM::UnknownQuantityHandler *RunawayFluid::GetUnknownHandlerSingleImpuritySpecies(DREAM::FVM::Grid *g, 
    const real_t IMPURITY_DENSITY, const len_t IMPURITY_Z0, const len_t IMPURITY_Z, 
    const real_t HYDROGEN_DENSITY, const real_t T_cold) {
    DREAM::FVM::UnknownQuantityHandler *uqh = new DREAM::FVM::UnknownQuantityHandler();

    len_t N_IONS = 2;
    len_t Z_IONS[2] = {1, IMPURITY_Z};

    len_t nZ0 = 0;
    for (len_t i = 0; i < N_IONS; i++)
        nZ0 += Z_IONS[i] + 1;

    this->id_ions = uqh->InsertUnknown(DREAM::OptionConstants::UQTY_ION_SPECIES, "0", g, nZ0);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_COLD,  "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_HOT,   "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_TOT,   "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_T_COLD,  "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_F_HOT,   "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_E_FIELD, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_J_TOT,   "0", g);

    real_t ni;
    // Set initial values
    const len_t N = nZ0*g->GetNr();
    real_t *nions = new real_t[N];
    real_t ncold = 0;
    real_t ntot = 0;
    len_t ionOffset = 0, rOffset = 0;
    for (len_t iIon = 0; iIon < N_IONS; iIon++) 
        for (len_t Z0 = 0; Z0 <= Z_IONS[iIon]; Z0++, ionOffset++){
            if ((iIon == 0) && (Z0 ==1 )) // I know this isn't the most clever way, but it required very little brain power... 
                ni = HYDROGEN_DENSITY; 
            else if ((iIon == 1) && (Z0 == IMPURITY_Z0))
                ni = IMPURITY_DENSITY;
            else
                ni = 0;
            
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
    vector<string> tritiumNames(0);
    vector<string> names(N_IONS);
    for (len_t i = 0; i < N_IONS; i++)
        names[i] = "";//ION_NAMES[i];

    return new DREAM::IonHandler(
        g->GetRadialGrid(), uqh, Z_IONS, N_IONS, names, tritiumNames
    );
}

DREAM::RunawayFluid *RunawayFluid::GetRunawayFluid(
    DREAM::CollisionQuantity::collqty_settings *cq, const len_t N_IONS,
    const len_t *Z_IONS, const real_t ION_DENSITY_REF, const real_t T_cold,
    const real_t B0, const len_t nr, 
    enum DREAM::OptionConstants::eqterm_dreicer_mode dreicer_mode,
    enum DREAM::OptionConstants::collqty_Eceff_mode eceff_mode
){
    DREAM::FVM::Grid *grid = this->InitializeFluidGrid(nr,B0);
    DREAM::FVM::UnknownQuantityHandler *unknowns = GetUnknownHandler(grid,N_IONS, Z_IONS, ION_DENSITY_REF,T_cold);

    DREAM::IonHandler *ionHandler = GetIonHandler(grid,unknowns, N_IONS, Z_IONS);
    ionHandler->Rebuild();
    DREAM::OptionConstants::momentumgrid_type gridtype = DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI;

    DREAM::CoulombLogarithm *lnLEE = new DREAM::CoulombLogarithm(grid,unknowns,ionHandler,gridtype,cq,DREAM::CollisionQuantity::LNLAMBDATYPE_EE);
    DREAM::CoulombLogarithm *lnLEI = new DREAM::CoulombLogarithm(grid,unknowns,ionHandler,gridtype,cq,DREAM::CollisionQuantity::LNLAMBDATYPE_EI);
    DREAM::SlowingDownFrequency *nuS = new DREAM::SlowingDownFrequency(grid,unknowns,ionHandler,lnLEE,lnLEI,gridtype,cq);
    DREAM::PitchScatterFrequency *nuD = new DREAM::PitchScatterFrequency(grid,unknowns,ionHandler,lnLEI,lnLEE,gridtype,cq);

    DREAM::RunawayFluid *REFluid = new DREAM::RunawayFluid(grid, unknowns, nuS,nuD,lnLEE,lnLEI, cq, ionHandler, DREAM::OptionConstants::CONDUCTIVITY_MODE_BRAAMS, dreicer_mode, eceff_mode, DREAM::OptionConstants::EQTERM_AVALANCHE_MODE_FLUID, DREAM::OptionConstants::EQTERM_COMPTON_MODE_NEGLECT, 0.0);
    REFluid->Rebuild();
    return REFluid;
}

DREAM::RunawayFluid *RunawayFluid::GetRunawayFluidSingleImpuritySpecies(
    DREAM::CollisionQuantity::collqty_settings *cq, const real_t IMPURITY_DENSITY,
    const len_t IMPURITY_Z0, const len_t IMPURITY_Z,
    const real_t B0, 
    enum DREAM::OptionConstants::eqterm_dreicer_mode dreicer_mode,
    enum DREAM::OptionConstants::collqty_Eceff_mode eceff_mode,
    const real_t HYDROGEN_DENSITY, const real_t T_cold
){
    len_t nr  = 1;
    DREAM::FVM::Grid *grid = this->InitializeFluidGrid(nr,B0);
    DREAM::FVM::UnknownQuantityHandler *unknowns = GetUnknownHandlerSingleImpuritySpecies(grid, IMPURITY_DENSITY, IMPURITY_Z0, IMPURITY_Z,
       HYDROGEN_DENSITY, T_cold);
    
    len_t N_IONS = 2; len_t Z_IONS[2] = {1, IMPURITY_Z};
    DREAM::IonHandler *ionHandler = GetIonHandler(grid,unknowns, N_IONS, Z_IONS);
    ionHandler->Rebuild();
    DREAM::OptionConstants::momentumgrid_type gridtype = DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI;

    DREAM::CoulombLogarithm *lnLEE = new DREAM::CoulombLogarithm(grid,unknowns,ionHandler,gridtype,cq,DREAM::CollisionQuantity::LNLAMBDATYPE_EE);
    DREAM::CoulombLogarithm *lnLEI = new DREAM::CoulombLogarithm(grid,unknowns,ionHandler,gridtype,cq,DREAM::CollisionQuantity::LNLAMBDATYPE_EI);
    DREAM::SlowingDownFrequency *nuS = new DREAM::SlowingDownFrequency(grid,unknowns,ionHandler,lnLEE,lnLEI,gridtype,cq);
    DREAM::PitchScatterFrequency *nuD = new DREAM::PitchScatterFrequency(grid,unknowns,ionHandler,lnLEI,lnLEE,gridtype,cq);

    DREAM::RunawayFluid *REFluid = new DREAM::RunawayFluid(grid, unknowns, nuS,nuD,lnLEE,lnLEI, cq, ionHandler, DREAM::OptionConstants::CONDUCTIVITY_MODE_BRAAMS, dreicer_mode, eceff_mode, DREAM::OptionConstants::EQTERM_AVALANCHE_MODE_FLUID, DREAM::OptionConstants::EQTERM_COMPTON_MODE_NEGLECT, 0.0);
    REFluid->Rebuild();
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

    // This first test (until line 321) could maybe be removed now, but I guess more testing doesn't do any harm
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
    T_cold = 50.0; // eV
    B0 = 3.0;
    REFluid = GetRunawayFluid(cq,N_IONS2, Z_IONS2, ION_DENSITY_REF, T_cold,B0,nr);
    Eceff3 = REFluid->GetEffectiveCriticalField(0);

    real_t TabulatedEceff1 = 8.88081;
    real_t TabulatedEceff2 = 8.00666;
    real_t TabulatedEceff3 = 1.10307;
    real_t delta1 = abs(Eceff1-TabulatedEceff1)/TabulatedEceff1;
    real_t delta2 = abs(Eceff2-TabulatedEceff2)/TabulatedEceff2;
    real_t delta3 = abs(Eceff3-TabulatedEceff3)/TabulatedEceff3;
    real_t threshold = 1e-3; 
    bool success = (delta1 < threshold) && (delta2 < threshold) && (delta3 < threshold);

    // Next part of the test, used to target the PPCF implementation. The plasma composition is chosen from the paper, but compared with numerical values from the script on GitHub (with He)
    constexpr int_t N_PLASMAS_TO_TEST = 5;
    constexpr int_t N_MODES = 3;
    DREAM::OptionConstants::eqterm_dreicer_mode dm=DREAM::OptionConstants::EQTERM_DREICER_MODE_NONE; // need to set dreicer in order to set Eceff mode
    DREAM::OptionConstants::collqty_Eceff_mode ECEFF_MODES[N_MODES] =  
        {DREAM::OptionConstants::COLLQTY_ECEFF_MODE_CYLINDRICAL, DREAM::OptionConstants::COLLQTY_ECEFF_MODE_FULL, DREAM::OptionConstants::COLLQTY_ECEFF_MODE_SIMPLE};
    
    len_t  Z_IMPURITY[N_PLASMAS_TO_TEST]               = { 18,    18,   18,   10,    2};
    len_t  Z0_IMPURITY[N_PLASMAS_TO_TEST]              = {  1,     1,    4,    1,    2};
    real_t B0_LIST[N_PLASMAS_TO_TEST]                  = {0.1,     5,  0.1,  0.1,    5};
    real_t IMPURITY_DENSITY[N_PLASMAS_TO_TEST]         = {1e20, 1e20, 1e19, 1e20, 1e21};
    real_t Eceff;
    real_t ECEFF_TABULATED_2[N_MODES][N_PLASMAS_TO_TEST] = {{1.75462, 2.04106, 0.27224, 0.88817, 2.14834},
                                                         //  1.65448, 1.97123, 0.25948, 0.85776, 2.10482, N_A = 100 sampling 1/x -1
                                                            {1.65449, 1.97124, 0.25948, 0.85776, 2.10482}, // prev implementation (without splines)
                                                            {1.65449, 1.97124, 0.25948, 0.85776, 2.10482} // TODO: update
                                                            }; 

    real_t delta; 
    for (len_t eceffMode = 0; eceffMode<N_MODES; eceffMode++){
        for (len_t i_test=0; i_test< N_PLASMAS_TO_TEST; i_test++){
            REFluid = GetRunawayFluidSingleImpuritySpecies(cq, IMPURITY_DENSITY[i_test],
                Z0_IMPURITY[i_test], Z_IMPURITY[i_test], B0_LIST[i_test], dm, ECEFF_MODES[eceffMode]);
            Eceff = REFluid->GetEffectiveCriticalField(0);
            delta = abs(Eceff-ECEFF_TABULATED_2[eceffMode][i_test])/ECEFF_TABULATED_2[eceffMode][i_test];
            success = success && delta < threshold;
        }
    }
    delete REFluid;

    return success;
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

    #define NR 3
    len_t nr = NR;
    const len_t N_IONS = 2;
    const len_t Z_IONS[N_IONS] = {10,18};
    real_t ION_DENSITY_REF = 1e18; // m-3
    real_t T_cold = 1; // eV
    real_t B0 = 5;
    DREAM::RunawayFluid *REFluid = GetRunawayFluid(cq,N_IONS, Z_IONS, ION_DENSITY_REF, T_cold,B0,nr);

    const real_t *GammaAva =  REFluid->GetAvalancheGrowthRate();
    const real_t GammaTabulated[NR] = {159.791, 11533.2, 24326.7};

    real_t *deltas = new real_t[NR];
    for(len_t ir=0; ir<nr;ir++)
        deltas[ir] = abs(GammaAva[ir]-GammaTabulated[ir])/GammaTabulated[ir];

/*
    cout << "Delta1: " << deltas[1] << endl;
    cout << "Delta2: " << deltas[2] << endl;
    cout << "Delta3: " << deltas[3] << endl;
*/
    real_t threshold = 2e-2;

    bool success = true;
    for (len_t ir = 0; ir < nr; ir++) {
        if (deltas[ir] > threshold) {
            this->PrintError(
                "Avalanche growth-rate deviates from tabulated values at ir = "
                LEN_T_PRINTF_FMT ".",
                ir
            );
            success = false;
            break;
        }
    }

    delete [] deltas;
    delete cq;
    delete REFluid;

    return success;

    #undef NR
}

real_t RunawayFluid::_ConnorHastieFormula(
    const real_t ne, const real_t EED,
    const real_t EEc, const real_t Zeff, const real_t tauEE,
    bool withCorrections
) {
    /* GO implementation
     * -----------------
    Te = electron.T;
    ne = electron.n;
    
    u2=Te*PhysConst.qe/(PhysConst.me*PhysConst.vlight^2) ./ (ne./(discr.nei0));
    
    %Correct the loglambda in E_C used to normalize the
    %electric field (to get the correct E/E_D)
    EoED = EfieldA.*u2./(electron.loglambda./tokamak.lnLambda);
    
    EoverEC=EfieldA*discr.nei0./electron.n.*tokamak.lnLambda./electron.loglambda_c;
    if(eta_lambda_h_correction)
        lambda=8*EoverEC.*(EoverEC-1/2-sqrt(EoverEC.^2-EoverEC));
        eta=1/4*EoverEC.^2./(EoverEC-1).*(pi/2-asin(1-2./EoverEC)).^2;
        h=1./(3*((EoverEC-1).^3).^1/2).*(((EoverEC-1).^3).^1/2.*8./(Zeff+1)+sqrt(EoverEC-1).*(1+2*sqrt(EoverEC))-2);
        
        eta(EoverEC<1)=0;%EoverEC<1 causes eta to give a complex value
        %which causes problems. F(EoverEC<1)=0 anyway, so we can might as
        %well set eta to zero in these cases.
    else
        lambda=1;
        eta=1;
        h=1;
    end
    F=3*0.99*electron.loglambda/2/sqrt(pi)*discr.nei0*PhysConst.qe*PhysConst.vlight...
        /discr.Ji0.*sqrt(ne/discr.nei0).*...
        u2.^(-3/2).*EoED.^(-3/16*(1+Zeff).*h).*...
        exp( -lambda.*1/4./EoED - sqrt(eta.*(1+Zeff)./EoED) );
    
    F(EoverEC<1)=0;%No runaways can be generated when E<EC.
    */

    real_t lambda=1, eta=1, h=1;

    if (withCorrections) {
        lambda = 8*EEc*(EEc-1.0/2.0-sqrt(EEc*EEc-EEc));
        eta    = 1.0/4.0*EEc*EEc/(EEc-1)*
                    (M_PI/2.0-asin(1-2./EEc))*
                    (M_PI/2.0-asin(1-2./EEc));
        h=1./(3*sqrt(EEc-1)*(EEc-1))*
            (sqrt(EEc-1)*(EEc-1)*8./(Zeff+1)+sqrt(EEc-1)*(1+2*sqrt(EEc))-2);

        this->PrintStatus("[TS] lambda = %.16e", lambda);
        this->PrintStatus("[TS] eta    = %.16e", eta);
        this->PrintStatus("[TS] h      = %.16e", h);
    }

    real_t C=1;
    return 
        C*ne/tauEE*
        pow(EED, (-3.0/16.0*(1+Zeff)*h))*
        exp( -lambda/(4.*EED) - sqrt(eta*(1+Zeff)/EED) );
}

/**
 * Compare the Connor-Hastie Dreicer runaway rate with tabulated values.
 */
bool RunawayFluid::CompareConnorHastieRateWithTabulated() {
    bool success = true;
    DREAM::CollisionQuantity::collqty_settings *cq =
        new DREAM::CollisionQuantity::collqty_settings;

    cq->collfreq_type = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED;
    cq->collfreq_mode = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
    cq->lnL_type      = DREAM::OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;
    cq->bremsstrahlung_mode = DREAM::OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT;
    cq->pstar_mode = DREAM::OptionConstants::COLLQTY_PSTAR_MODE_COLLISIONLESS;

    #define NR 3
    len_t nr = NR;
    const len_t N_IONS = 2;
    const len_t Z_IONS[N_IONS] = {10,18};
    real_t ION_DENSITY_REF = 1e18; // m-3
    real_t T_cold = 300; // eV
    real_t B0 = 5;
    DREAM::RunawayFluid *REFluid = GetRunawayFluid(cq,N_IONS, Z_IONS, ION_DENSITY_REF, T_cold,B0,nr, DREAM::OptionConstants::EQTERM_DREICER_MODE_CONNOR_HASTIE_NOCORR);

    DREAM::FVM::UnknownQuantityHandler *uqn = REFluid->GetUnknowns();
    len_t id_n_cold = uqn->GetUnknownID(DREAM::OptionConstants::UQTY_N_COLD);
    real_t ncold = uqn->GetUnknownData(id_n_cold)[0];
    this->PrintStatus("ncold = %e", ncold);
    //real_t ncold = 4.6009999999999997e+21;
    real_t Zeff  = REFluid->GetIonHandler()->GetZeff(0);

    real_t Ec = REFluid->GetConnorHastieField_COMPLETESCREENING(0);
    real_t ED = REFluid->GetDreicerElectricField(0);
    real_t tauEE = REFluid->GetElectronCollisionTimeThermal(0);

    real_t Emin = 0.001, Emax = 0.1;
    len_t nE = 10;
    bool withCorrections = false;

    // Evaluate Connor-Hastie rate
    DREAM::ConnorHastie *ch = REFluid->GetConnorHastieRunawayRate();
    ch->IncludeCorrections(withCorrections);

    const real_t TOLERANCE = 100.0*std::numeric_limits<real_t>::epsilon();
    for (len_t i = 0; i < nE; i++) {
        //real_t E = 0.02*ED;
        real_t E = Emin + i*(Emax-Emin)/(nE-1);

        real_t DREAMrate = ch->RunawayRate(0, E, ncold, Zeff);
        real_t GOrate    = _ConnorHastieFormula(ncold, E/ED, E/Ec, Zeff, tauEE, withCorrections);

        real_t Delta;
        if (GOrate == 0)
            Delta = abs(DREAMrate);
        else
            Delta = abs((DREAMrate-GOrate)/GOrate);

        if (Delta > TOLERANCE) {
            this->PrintError(
                "DREAM and GO Connor-Hastie runaway rates do not agree at E = %e. Delta = %e",
                E, Delta
            );
            success = false;
            break;
        }
    }

    delete cq;
    delete REFluid;
    delete uqn;

    return success;

    #undef NR
}

