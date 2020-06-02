
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

#include "gsl/gsl_integration.h"

using namespace DREAMTESTS::FVM;
using namespace std;






/**
 * Run this test.
 */
bool AnalyticBRadialGridGenerator::Run(bool) {
    Initialize();

    bool success = true;
    if (TestGeneralFluxSurfaceAverage())
        this->PrintOK("The general flux surface average agrees with calculated reference value.");
    else {
        success = false;
        this->PrintError("General flux surface average test failed.");
    }

    if (TestGeneralBounceAverage())
        this->PrintOK("The general bounce average agrees with calculated reference value.");
    else {
        success = false;
        this->PrintError("General bounce average test failed.");
    }

    if (CompareBounceAverageMethods())
        this->PrintOK("The 2 different methods of evaluating the general bounce average agree with eachother.");
    else {
        success = false;
        this->PrintError("General bounce average test failed.");
    }

    

    delete grid;

    return success;
}


    /**
     * Generates arbitrary profiles that at ir=1 corresponds to:
     * r          = 1.1186499999999999
     * G          = 3.3403978996084018
     * GPrime     = 2.818765586918027
     * kappa      = 2.4409505
     * kappaPrime = 1.696103135079847
     * delta      = 0.22655
     * deltaPrime = 0.61975106168968697
     * Delta      = 0.23591500000000001
     * DeltaPrime = 0.64536999200088296
     * psi        = 0.69564999999999999
     * psiPrime   = 1.9030228474860778
     */
void AnalyticBRadialGridGenerator::Initialize(){

    const len_t np = 20, nxi = 5;//nxi = 15;
    const real_t pMin = 0, pMax = 10;
    len_t nrProfiles = 53;
//    len_t ntheta_ref = 10000;
    len_t ntheta_ref = 5000;
//    len_t ntheta_interp = 100;
    len_t ntheta_interp = 500;
    len_t nr = 3;

    real_t r0 = 0.7531;
    real_t ra = 1.4842;
    real_t R0 = 3.39431;
    real_t 
        *rProfiles = new real_t[nrProfiles],
        *Gs = new real_t[nrProfiles], 
        *psi_p0s = new real_t[nrProfiles], 
        *kappas = new real_t[nrProfiles], 
        *deltas = new real_t[nrProfiles], 
        *Deltas = new real_t[nrProfiles];

    for (len_t it = 0; it<nrProfiles; it++){
        rProfiles[it] = r0 + it*(ra-r0)/(nrProfiles-1);
        Gs[it]      = (2.8353 + 2.1828*it*(it-1)/(nrProfiles*nrProfiles))/R0; 
        kappas[it]  = 1.82094 + 1.240021*it/(nrProfiles-1);     
        Deltas[it]  = 0.47183*it/(nrProfiles-1);
        psi_p0s[it] = 1.3913*it/(nrProfiles-1);        
//        deltas[it]  = 1.9531*it/(nrProfiles-1);
        deltas[it]  = 0.4531*it/(nrProfiles-1);
    }

    auto *ABrgg = new DREAM::FVM::AnalyticBRadialGridGenerator(nr, r0, ra, R0, 
                        ntheta_ref, ntheta_interp, rProfiles, nrProfiles, Gs, psi_p0s,
                        kappas, deltas, Deltas); 


    auto *rg   = new DREAM::FVM::RadialGrid(ABrgg);


    auto *pgg = new DREAM::FVM::PXiGrid::PUniformGridGenerator(np, pMin, pMax);
    auto *xgg = new DREAM::FVM::PXiGrid::XiUniformGridGenerator(nxi);

    auto *mgg = new DREAM::FVM::PXiGrid::MomentumGridGenerator(pgg, xgg);
    auto *mg  = new DREAM::FVM::PXiGrid::PXiMomentumGrid(mgg, 0, rg);

    grid = new DREAM::FVM::Grid(rg, mg);

    grid->RebuildJacobians();
    
}

bool AnalyticBRadialGridGenerator::CompareBounceAverageMethods(){
    bool success = false;

    std::function<real_t(real_t,real_t)> 
        generalFunction = [](real_t x, real_t y)
        {return 0.315*pow(abs(x),1.382)*pow(y,2.913);} ;
    
    len_t ir = 1;
    len_t i = 10;
    len_t j = 3;
//    len_t j = 14;

    DREAM::FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
    real_t generalBounceAverage = grid->GetRadialGrid()->CalculateBounceAverage(mg,ir,i,j,DREAM::FVM::FLUXGRIDTYPE_P1,generalFunction);
//    real_t r = grid->GetRadialGrid()->GetR(ir);
    real_t p = grid->GetMomentumGrid(ir)->GetP_f1(i,j);
    real_t xi0 = grid->GetMomentumGrid(ir)->GetXi0_f1(i,j);    
    gsl_integration_workspace *gsl_ad_w = gsl_integration_workspace_alloc(1000);
    real_t bounceAverageAtP = grid->GetRadialGrid()->evaluatePXiBounceAverageAtP(ir,p,xi0,false,generalFunction,gsl_ad_w);
    real_t relativeError = abs(bounceAverageAtP-generalBounceAverage)/bounceAverageAtP;

    /*
    cout << "r: " << r << ", p: " << p << ", xi0: " << xi0 << endl;
    cout << "CalculateBounceAverage: " << generalBounceAverage << endl;
    cout << "evaluateBounceAverageAtP: " << bounceAverageAtP << endl;
    cout << "Relative error: " << relativeError << endl;
    */
    if (relativeError < 1e-3)
        success = true;

    return success;

}


bool AnalyticBRadialGridGenerator::TestGeneralBounceAverage(){
    bool success = false;

    std::function<real_t(real_t,real_t)> 
        generalFunction = [](real_t x, real_t y)
        {return 0.315*pow(abs(x),1.382)*pow(y,2.913);} ;
    
    len_t ir = 1;
    len_t i = 10;
    len_t j = 3;
//    len_t j = 10;

    DREAM::FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
    real_t generalBounceAverage = grid->GetRadialGrid()->CalculateBounceAverage(mg,ir,i,j,DREAM::FVM::FLUXGRIDTYPE_P1,generalFunction);
    real_t referenceValueMatlab = 0.212990959899721;
    
    real_t relativeError = abs(referenceValueMatlab-generalBounceAverage)/referenceValueMatlab;
    //cout << "Relative error: " << relativeError << endl;

    if (relativeError < 1e-3)
        success = true;

    return success;
}

bool AnalyticBRadialGridGenerator::TestGeneralFluxSurfaceAverage(){
    bool success = false;

    std::function<real_t(real_t,real_t,real_t)> 
        generalFunction = [](real_t x, real_t y, real_t z)
        {return .13153*pow(x,1.4313)*pow(y,0.3901)*pow(z,2.159);} ;
    real_t generalFluxSurfaceAverage = grid->GetRadialGrid()->CalculateFluxSurfaceAverage(1,false,generalFunction);
    //real_t referenceValueMatlab = 2.738968863114242;
    real_t referenceValueMatlab = 0.791833837394785;

    real_t relativeError =  abs(generalFluxSurfaceAverage - referenceValueMatlab)/referenceValueMatlab;
    //cout << "Relative error: " << relativeError << endl; 
    
    if (relativeError < 1e-5)
        success = true;
    return success;
}
