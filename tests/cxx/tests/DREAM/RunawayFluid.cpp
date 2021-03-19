/**
 * A test for RunawayFluid in DREAM, that verifies the implementation of 
 * derived collision quantities such as the avalanche growth rate and critical E field. 
 * Calculations are benchmarked with values tabulated by DREAM simulations in
 * commit c8f1923d962b3b565ace4e2b033e37ad0a0cb5a8.
 * The Eceff values were updated in commit f69903e9aa0d469b55dfd8e6f06e27343b871e50
 * when a minor bug in the bremsstrahlung force was fixed which caused a <1% error.
 *  
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

    if (VerifyAnalyticalDistributionRE())
        this->PrintOK("The analytical RE pitch distribution passes all tests.");
    else {
        success = false;
        this->PrintError("The analytical RE pitch distribution test failed.");
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
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_COLD,  "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_HOT,   "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_RE,    "0", g);
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
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_RE,    "0", g);
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
    const len_t N_IONS, const len_t *Z_IONS, const real_t ION_DENSITY_REF, 
    const real_t T_cold, const real_t B0, const len_t nr, bool generalGrid,
    enum DREAM::OptionConstants::eqterm_dreicer_mode dreicer_mode,
    enum DREAM::OptionConstants::collqty_Eceff_mode eceff_mode
){
    DREAM::FVM::Grid *grid;
    if(generalGrid)
        grid = this->InitializeGridGeneralFluid(nr);
    else 
        grid = this->InitializeFluidGrid(nr,B0);
    
    DREAM::FVM::UnknownQuantityHandler *unknowns = GetUnknownHandler(grid,N_IONS, Z_IONS, ION_DENSITY_REF,T_cold);

    DREAM::IonHandler *ionHandler = GetIonHandler(grid,unknowns, N_IONS, Z_IONS);
    ionHandler->Rebuild();

    return ConstructRunawayFluid(grid, unknowns, ionHandler,dreicer_mode, eceff_mode);
}

DREAM::RunawayFluid *RunawayFluid::GetRunawayFluidSingleImpuritySpecies(
    const real_t IMPURITY_DENSITY,
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

    return ConstructRunawayFluid(grid, unknowns, ionHandler,dreicer_mode, eceff_mode);
}

DREAM::RunawayFluid *RunawayFluid::ConstructRunawayFluid(
    DREAM::FVM::Grid *grid, DREAM::FVM::UnknownQuantityHandler *unknowns, 
    DREAM::IonHandler *ionHandler, enum DREAM::OptionConstants::eqterm_dreicer_mode dreicer_mode,
    enum DREAM::OptionConstants::collqty_Eceff_mode eceff_mode
) {
    DREAM::CollisionQuantity::collqty_settings
        *cqPc = new DREAM::CollisionQuantity::collqty_settings,
        *cqEc = new DREAM::CollisionQuantity::collqty_settings;
    cqPc->collfreq_type = cqEc->collfreq_type = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED;
    cqPc->collfreq_mode = cqEc->collfreq_mode = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
    cqPc->lnL_type      = cqEc->lnL_type      = DREAM::OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;
    cqPc->pstar_mode    = cqEc->pstar_mode    = DREAM::OptionConstants::COLLQTY_PSTAR_MODE_COLLISIONLESS;
    cqPc->bremsstrahlung_mode = DREAM::OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT;
    cqEc->bremsstrahlung_mode = DREAM::OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_STOPPING_POWER;

    DREAM::OptionConstants::momentumgrid_type gridtype = DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI;

    DREAM::CoulombLogarithm *lnLEE = new DREAM::CoulombLogarithm(grid,unknowns,ionHandler,gridtype,cqPc,DREAM::CollisionQuantity::LNLAMBDATYPE_EE);
    DREAM::CoulombLogarithm *lnLEI = new DREAM::CoulombLogarithm(grid,unknowns,ionHandler,gridtype,cqPc,DREAM::CollisionQuantity::LNLAMBDATYPE_EI);
    DREAM::SlowingDownFrequency *nuS  = new DREAM::SlowingDownFrequency(grid,unknowns,ionHandler,lnLEE,lnLEI,gridtype,cqPc);
    DREAM::PitchScatterFrequency *nuD = new DREAM::PitchScatterFrequency(grid,unknowns,ionHandler,lnLEI,lnLEE,gridtype,cqPc);

    DREAM::AnalyticDistributionRE::dist_mode re_dist_mode = (eceff_mode==DREAM::OptionConstants::COLLQTY_ECEFF_MODE_SIMPLE) ? 
            DREAM::AnalyticDistributionRE::RE_PITCH_DIST_SIMPLE : DREAM::AnalyticDistributionRE::RE_PITCH_DIST_FULL;
    DREAM::AnalyticDistributionRE *distRE =  new DREAM::AnalyticDistributionRE(grid->GetRadialGrid(), unknowns, nuD, cqEc, re_dist_mode, 100*sqrt(std::numeric_limits<real_t>::epsilon()));
    DREAM::RunawayFluid *REFluid = new DREAM::RunawayFluid(
        grid, unknowns, nuS, nuD,lnLEE,lnLEI, ionHandler, distRE, cqPc, cqEc,
        DREAM::OptionConstants::CONDUCTIVITY_MODE_BRAAMS, dreicer_mode, 
        eceff_mode, DREAM::OptionConstants::EQTERM_AVALANCHE_MODE_FLUID, 
        DREAM::OptionConstants::EQTERM_COMPTON_MODE_NEGLECT, 0.0
    );
    REFluid->Rebuild();
    return REFluid;
}

bool RunawayFluid::CompareEceffWithTabulated(){
    len_t nr = 1;

    // This first test (until line 321) could maybe be removed now, but I guess more testing doesn't do any harm
    real_t Eceff1, Eceff2, Eceff3;
    const len_t N_IONS = 2;
    const len_t Z_IONS[N_IONS] = {10,18};
    real_t ION_DENSITY_REF = 1e18; // m-3
    real_t T_cold = 1; // eV
    real_t B0 = 5;

    DREAM::RunawayFluid *REFluid = GetRunawayFluid(N_IONS, Z_IONS, ION_DENSITY_REF, T_cold,B0,nr);
    Eceff1 = REFluid->GetEffectiveCriticalField(0);
    delete REFluid;

    B0 = 0.1;
    REFluid = GetRunawayFluid(N_IONS, Z_IONS, ION_DENSITY_REF, T_cold,B0,nr);
    Eceff2 = REFluid->GetEffectiveCriticalField(0);
    delete REFluid;

    const len_t N_IONS2 = 1;
    const len_t Z_IONS2[N_IONS2] = {2};
    ION_DENSITY_REF = 1e20; // m-3
    T_cold = 50.0; // eV
    B0 = 3.0;
    REFluid = GetRunawayFluid(N_IONS2, Z_IONS2, ION_DENSITY_REF, T_cold,B0,nr);
    Eceff3 = REFluid->GetEffectiveCriticalField(0);
    delete REFluid;

    real_t TabulatedEceff1 = 8.928710;
    real_t TabulatedEceff2 = 8.062882;
    real_t TabulatedEceff3 = 1.103813;
    real_t delta1 = abs(Eceff1-TabulatedEceff1)/TabulatedEceff1;
    real_t delta2 = abs(Eceff2-TabulatedEceff2)/TabulatedEceff2;
    real_t delta3 = abs(Eceff3-TabulatedEceff3)/TabulatedEceff3;
    real_t threshold = 1e-3; 
    bool success = (delta1 < threshold) && (delta2 < threshold) && (delta3 < threshold);
    if(delta1>threshold)
        this->PrintWarning("Eceff Delta 1: %e", delta1);
    if(delta2>threshold)
        this->PrintWarning("Eceff Delta 2: %e", delta2);
    if(delta3>threshold)
        this->PrintWarning("Eceff Delta 3: %e", delta3);
     
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
    real_t ECEFF_TABULATED_2[N_MODES][N_PLASMAS_TO_TEST] = {{1.754624, 2.041058, 0.272244, 0.888174, 2.148342},
                                                            {1.666536, 1.980693, 0.260755, 0.861354, 2.106067},
                                                            {1.666536, 1.980693, 0.260755, 0.861354, 2.106067} 
                                                            }; 

    real_t delta; 
    for (len_t eceffMode = 0; eceffMode<N_MODES; eceffMode++){
        for (len_t i_test=0; i_test< N_PLASMAS_TO_TEST; i_test++){
            REFluid = GetRunawayFluidSingleImpuritySpecies(IMPURITY_DENSITY[i_test],
                Z0_IMPURITY[i_test], Z_IMPURITY[i_test], B0_LIST[i_test], dm, ECEFF_MODES[eceffMode]);
            Eceff = REFluid->GetEffectiveCriticalField(0);
            delta = abs(Eceff-ECEFF_TABULATED_2[eceffMode][i_test])/ECEFF_TABULATED_2[eceffMode][i_test];
            if(delta>threshold)
                this->PrintWarning("Eceff mode test Delta: %e", delta);
            success = success && delta < threshold;
            delete REFluid;
        }
    }
    return success;
}

/**
 * Evalutes the semi-analytic avalanche growth rate in a Neon-Argon plasma for 
 * three different E fields and compares with tabulated values.
 */
bool RunawayFluid::CompareGammaAvaWithTabulated(){
    #define NR 3
    len_t nr = NR;
    const len_t N_IONS = 2;
    const len_t Z_IONS[N_IONS] = {10,18};
    real_t ION_DENSITY_REF = 1e18; // m-3
    real_t T_cold = 1; // eV
    real_t B0 = 5;
    DREAM::RunawayFluid *REFluid = GetRunawayFluid(N_IONS, Z_IONS, ION_DENSITY_REF, T_cold,B0,nr);

    const real_t *GammaAva =  REFluid->GetAvalancheGrowthRate();
    const real_t GammaTabulated[NR] = {159.791, 11533.2, 24326.7};

    real_t *deltas = new real_t[NR];
    for(len_t ir=0; ir<nr;ir++)
        deltas[ir] = abs(GammaAva[ir]-GammaTabulated[ir])/GammaTabulated[ir];

    real_t threshold = 2e-2;

    bool success = true;
    for (len_t ir = 0; ir < nr; ir++) {
        if (deltas[ir] > threshold) {
            this->PrintError(
                "Avalanche growth-rate deviates from tabulated values at ir = "
                LEN_T_PRINTF_FMT ". Delta = %f and GammaAva = %f",
                ir, deltas[ir], GammaAva[ir]
            );
            success = false;
            break;
        }
    }

    delete [] deltas;
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

    #define NR 3
    len_t nr = NR;
    const len_t N_IONS = 2;
    const len_t Z_IONS[N_IONS] = {10,18};
    real_t ION_DENSITY_REF = 1e18; // m-3
    real_t T_cold = 300; // eV
    real_t B0 = 5;
    DREAM::RunawayFluid *REFluid = GetRunawayFluid(N_IONS, Z_IONS, ION_DENSITY_REF, T_cold,B0,nr, false, DREAM::OptionConstants::EQTERM_DREICER_MODE_CONNOR_HASTIE_NOCORR);

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

    delete REFluid;
    delete uqn;

    return success;

    #undef NR
}


/**
 * Verify tractable analytical limits of the RE distribution-averaged equation terms  
 */
bool RunawayFluid::VerifyAnalyticalDistributionRE(){
    bool success = true;

    #define NR 3
    len_t nr = NR;
    const len_t N_IONS = 2;
    const len_t Z_IONS[N_IONS] = {10,18};
    real_t ION_DENSITY_REF = 1e18; // m-3
    real_t T_cold = 1; // eV
    real_t B0 = 5;
    DREAM::RunawayFluid *REFluid = GetRunawayFluid(N_IONS, Z_IONS, ION_DENSITY_REF, T_cold,B0,nr,true);
    DREAM::AnalyticDistributionRE *distRE = REFluid->GetAnalyticDistributionRE();
    DREAM::FVM::RadialGrid *rGrid = distRE->GetRadialGrid();

    // Generate distribution-averaged {A^p}-coefficients of the synchrotron and electric field equation terms 
    real_t synchrotronPrefactor = DREAM::Constants::ec * DREAM::Constants::ec * DREAM::Constants::ec * DREAM::Constants::ec 
            / ( 6.0 * M_PI * DREAM::Constants::eps0 * DREAM::Constants::me * DREAM::Constants::me * DREAM::Constants::me
            * DREAM::Constants::c * DREAM::Constants::c * DREAM::Constants::c);
    std::function<real_t(real_t)> synchPitchFunc = [](real_t xi0){return 1.0-xi0*xi0;}; 
    DREAM::REPitchDistributionAveragedBACoeff *AveragedSynchrotronTerm = new DREAM::REPitchDistributionAveragedBACoeff(
        rGrid, distRE, &(DREAM::FVM::RadialGrid::BA_FUNC_B_CUBED), nullptr, 
        DREAM::FVM::RadialGrid::BA_PARAM_B_CUBED, synchPitchFunc,
        [synchrotronPrefactor,rGrid](len_t ir, real_t p){
            return -p*sqrt(1+p*p)*synchrotronPrefactor*rGrid->GetBmin(ir)*rGrid->GetBmin(ir);
    });
    AveragedSynchrotronTerm->GridRebuilt();
    
    DREAM::REPitchDistributionAveragedBACoeff *AveragedEFieldTerm = new DREAM::REPitchDistributionAveragedBACoeff(
        rGrid, distRE, &(DREAM::FVM::RadialGrid::BA_FUNC_XI), nullptr, 
        DREAM::FVM::RadialGrid::BA_PARAM_XI,[](real_t xi0){return xi0;},
        [rGrid](len_t ir, real_t){
            return DREAM::Constants::ec  / (DREAM::Constants::me * DREAM::Constants::c) 
                * sqrt(rGrid->GetFSA_B2(ir)) / rGrid->GetFSA_B(ir);
    });
    AveragedEFieldTerm->GridRebuilt();

    len_t ir = 1;
    real_t p = 0.8;

    real_t infA = std::numeric_limits<real_t>::infinity();
    real_t zeroA = 0;
    
    real_t synchPrefact = -p*sqrt(1+p*p)*synchrotronPrefactor*rGrid->GetBmin(ir)*rGrid->GetBmin(ir);
    real_t synchAvg1 = AveragedSynchrotronTerm->EvaluateREPitchDistAverage(ir, p, &infA);
    real_t synchAvg1Alt = synchPitchFunc(1.0)*synchPrefact*rGrid->CalculatePXiBounceAverageAtP(
        ir,1.0, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, &(DREAM::FVM::RadialGrid::BA_FUNC_B_CUBED), nullptr, 
        DREAM::FVM::RadialGrid::BA_PARAM_B_CUBED
    );
    real_t synchAvg2 = AveragedSynchrotronTerm->EvaluateREPitchDistAverage(ir, p, &zeroA);
    real_t synchAvg2Alt = (2.0/3.0)*synchPrefact * rGrid->CalculateFluxSurfaceAverage(ir, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, 
        DREAM::FVM::RadialGrid::FSA_FUNC_B_SQUARED, nullptr, DREAM::FVM::RadialGrid::FSA_PARAM_B_SQUARED);

    real_t EPrefact = DREAM::Constants::ec  / (DREAM::Constants::me * DREAM::Constants::c) 
                * sqrt(rGrid->GetFSA_B2(ir)) / rGrid->GetFSA_B(ir);
    real_t EAvg1 = AveragedEFieldTerm->EvaluateREPitchDistAverage(ir, p, &infA);
    real_t EAvg1Alt = EPrefact; 
    real_t EAvg2 = AveragedEFieldTerm->EvaluateREPitchDistAverage(ir, p, &zeroA);

    real_t TOLERANCE = 1e-3;
    real_t delta1 = fabs(synchAvg1-synchAvg1Alt)/synchAvg2Alt;
    real_t delta2 = fabs(synchAvg2-synchAvg2Alt)/synchAvg2Alt;

    real_t delta3 = fabs(EAvg1-EAvg1Alt)/EAvg1Alt;
    real_t delta4 = fabs(EAvg2)/EAvg1Alt;
    
     
     if(delta1 > TOLERANCE){
        this->PrintError("Averaged synchrotron term does not agree with analytical "
            "formula at A=inf. Delta = %f \n", delta1);
        success = false;
     }
     if(delta2 > TOLERANCE){
        this->PrintError("Averaged synchrotron term does not agree with analytical "
            "formula at A=0. Delta = %f \n", delta2);
        success = false;
     }
     if(delta3 > TOLERANCE){
        this->PrintError("Averaged E field term does not agree with analytical "
            "formula at A=inf. Delta = %f \n", delta3);
        success = false;
     }
     if(delta4 > TOLERANCE){
        this->PrintError("Averaged E field term does not agree with analytical "
            "formula at A=0. Delta = %f \n", delta4);
        success = false;
     }

    return success;
}
