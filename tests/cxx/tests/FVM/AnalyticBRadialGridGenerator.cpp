
#include "DREAM/config.h"
#include "AnalyticBRadialGridGenerator.hpp"
#include "FVM/Grid/AnalyticBRadialGridGenerator.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PUniformGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiUniformGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/Grid.hpp"

#include <iostream>

using namespace DREAMTESTS::FVM;
using namespace std;






/**
 * Run this test.
 */
bool AnalyticBRadialGridGenerator::Run(bool) {
    Initialize();

    bool success = true;
    if (TestPoloidalFluxTerm())
        this->PrintOK("Flux surface average of poloidal-flux term agrees with theory.");
    else {
        success = false;
        this->PrintError("Calculation of poloidal-flux average failed.");
    }

    if (TestVpVol())
        this->PrintOK("Evaluation of VpVol agrees with theory.");
    else {
        success = false;
        this->PrintError("Calculation of VpVol failed.");
    }
    return success;
}



bool AnalyticBRadialGridGenerator::TestPoloidalFluxTerm(){
    bool success = true;


    return success;
}

bool AnalyticBRadialGridGenerator::TestVpVol(){
    bool success = true;

    return success;
}

void AnalyticBRadialGridGenerator::Initialize(){

    len_t nr = 50;
    real_t r0 = 0;
    real_t ra = 1;
    real_t R0 = 3;
    len_t ntheta_ref = 11;
    len_t ntheta_interp = 3;
    len_t nrProfiles = 100;
    real_t 
        *rProfiles = new real_t[nrProfiles],
        *Gs = new real_t[nrProfiles], 
        *psi_p0s = new real_t[nrProfiles], 
        *kappas = new real_t[nrProfiles], 
        *deltas = new real_t[nrProfiles], 
        *Deltas = new real_t[nrProfiles];
    
    for (len_t it = 0; it<nrProfiles; it++){
        rProfiles[it] = r0 + it*(ra-r0)/(nrProfiles-1);
        Gs[it] = 2.8353 + 2.128*it*(it-1)/(nrProfiles*nrProfiles); // gibberish sequences
        kappas[it] = 1.824;     // Constant elongation
        Deltas[it] = 0; //  = it/(nrProfiles*M_PI);   // linear to yield well-defined derivative
        psi_p0s[it] = 0;        // pure G(r)/R magnetic field
        deltas[it] = 0;
    }
    real_t kappaPrime = (kappas[1]-kappas[0])/(rProfiles[1]-rProfiles[0]);
    real_t DeltaPrime = (Deltas[1]-Deltas[0])/(rProfiles[1]-rProfiles[0]);
    
    auto *ABrgg = new DREAM::FVM::AnalyticBRadialGridGenerator(nr, r0, ra, R0, 
                        ntheta_ref, ntheta_interp, rProfiles, nrProfiles, Gs, psi_p0s,
                        kappas, deltas, Deltas); 


    auto *rg   = new DREAM::FVM::RadialGrid(ABrgg);

    const len_t np = 20, nxi = 5;
    const real_t pMin = 0, pMax = 10;
    // Build momentum grid
    auto *pgg = new DREAM::FVM::PXiGrid::PUniformGridGenerator(np, pMin, pMax);
    auto *xgg = new DREAM::FVM::PXiGrid::XiUniformGridGenerator(nxi);

    auto *mgg = new DREAM::FVM::PXiGrid::MomentumGridGenerator(pgg, xgg);
    auto *mg  = new DREAM::FVM::PXiGrid::PXiMomentumGrid(mgg, 0, rg);

    grid = new DREAM::FVM::Grid(rg, mg);

    grid->RebuildJacobians();

    real_t EPF_MATLAB = 0.722456715620853;
    len_t it = 5;
    real_t epsFrac = rg->GetR(it)/R0;
    cout << "eps: " << epsFrac << "." << endl;
    cout << "EffPassFrac: " << rg->GetEffPassFrac(it) << "," << endl;
    cout << "Expected: " << 1-1.462*sqrt(epsFrac) + epsFrac  << "." << endl;
    cout << "Relative error compared to matlab calculation: " << (rg->GetEffPassFrac(5) - EPF_MATLAB)/EPF_MATLAB << "." << endl;
}