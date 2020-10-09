
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
    const len_t np = 3, nxi = 5;
    len_t nrProfiles = 53;
    len_t ntheta_ref = 5000;
    len_t ntheta_interp = 20;
    len_t nr = 3;

    grid = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_ref, ntheta_interp, nrProfiles);

    silentMode = true;   
}


bool AnalyticBRadialGridGenerator::CompareBounceAverageMethods(){
    // TODO: include the new arguments, 
    // with R/R0 and |nabla r|^2 dependence 
    std::function<real_t(real_t,real_t,real_t,real_t)> 
        generalFunction = [](real_t x, real_t y,real_t,real_t)
        {return 0.315*pow(abs(x),1.382)*pow(y,2.913);} ;
    
    len_t ir = 1;
    len_t i = 2;
    len_t j = 3;

    real_t generalBounceAverage = grid->CalculateBounceAverage(ir,i,j,DREAM::FVM::FLUXGRIDTYPE_P1,generalFunction);
    real_t p = grid->GetMomentumGrid(ir)->GetP_f1(i,j);
    real_t xi0 = grid->GetMomentumGrid(ir)->GetXi0_f1(i,j);    
    real_t bounceAverageAtP = grid->GetRadialGrid()->CalculatePXiBounceAverageAtP(ir,p,xi0,DREAM::FVM::FLUXGRIDTYPE_P1,generalFunction);
    real_t relativeError = abs(bounceAverageAtP-generalBounceAverage)/bounceAverageAtP;

    if(!silentMode){
        real_t r = grid->GetRadialGrid()->GetR(ir);
        cout << "CompareBounceAverageMethods:" << endl;
        cout << "----------------------------" << endl;
        cout << "r: " << r << ", p: " << p << ", xi0: " << xi0 << endl;
        cout << "CalculateBounceAverage: " << generalBounceAverage << endl;
        cout << "evaluateBounceAverageAtP: " << bounceAverageAtP << endl;
        cout << "Relative error: " << relativeError << endl;
    }

    return (relativeError < 1e-3);
}


bool AnalyticBRadialGridGenerator::TestGeneralBounceAverage(){
    // TODO: include the new arguments, 
    // with R/R0 and |nabla r|^2 dependence 
    std::function<real_t(real_t,real_t,real_t,real_t)> 
        generalFunction = [](real_t x, real_t y,real_t,real_t)
        {return 0.315*pow(abs(x),1.382)*pow(y,2.913);} ;
    
    len_t ir = 1;
    len_t i = 2;
    len_t j = 3;

    real_t generalBounceAverage = grid->CalculateBounceAverage(ir,i,j,DREAM::FVM::FLUXGRIDTYPE_P1,generalFunction);
    real_t referenceValueMatlab = 0.212990959899721;
    
    real_t relativeError = abs(referenceValueMatlab-generalBounceAverage)/referenceValueMatlab;
    
    if(!silentMode){
        cout << "TestGeneralBounceAverage:" << endl;
        cout << "-------------------------" << endl;
        cout << "General bounce average: " << generalBounceAverage << endl;
        cout << "Matlab reference value: " << referenceValueMatlab << endl;
        cout << "Relative error: " << relativeError << endl;
    }
    return (relativeError < 1e-3);
}

bool AnalyticBRadialGridGenerator::TestGeneralFluxSurfaceAverage(){
    std::function<real_t(real_t,real_t,real_t)> 
        generalFunction = [](real_t x, real_t y, real_t z)
        {return .13153*pow(x,1.4313)*pow(y,0.3901)*pow(z,2.159);} ;
    real_t generalFluxSurfaceAverage = grid->GetRadialGrid()->CalculateFluxSurfaceAverage(1,DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION,generalFunction);
    real_t referenceValueMatlab = 0.791833837394785;

    real_t relativeError =  abs(generalFluxSurfaceAverage - referenceValueMatlab)/referenceValueMatlab;

    if(!silentMode){
        cout << "TestFluxSurfaceAverage:" << endl;
        cout << "-------------------------" << endl;
        cout << "General flux surface average: " << generalFluxSurfaceAverage << endl;
        cout << "Matlab reference value: " << referenceValueMatlab << endl;
        cout << "Relative error: " << relativeError << endl; 
    }

    return (relativeError < 1e-5);
}
