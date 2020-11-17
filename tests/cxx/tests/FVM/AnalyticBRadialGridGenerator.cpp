
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
    len_t ntheta_interp = 30;
    len_t nr = 3;
    real_t pmin=0, pmax=10;
    grid = InitializeGridGeneralRPXi(
        nr, np, nxi, ntheta_interp, nrProfiles, pmin, pmax,
        DREAM::FVM::FluxSurfaceAverager::QUAD_FIXED_LEGENDRE,
        DREAM::FVM::FluxSurfaceAverager::QUAD_FIXED_CHEBYSHEV
    );
    grid_adaptive = InitializeGridGeneralRPXi(
        nr, np, nxi, ntheta_interp, nrProfiles, pmin, pmax,
        DREAM::FVM::FluxSurfaceAverager::QUAD_ADAPTIVE,
        DREAM::FVM::FluxSurfaceAverager::QUAD_ADAPTIVE
    );

    silentMode = true;   
}


real_t getRelativeError(real_t I1,real_t I2){
    real_t relativeError = 0;
    if(I1||I2)
        relativeError = fabs(I2-I1)/(I1+I2);
    return relativeError;
}

bool AnalyticBRadialGridGenerator::CompareBounceAverageMethods(){
    // TODO: include the new arguments, 
    // with R/R0 and |nabla r|^2 dependence 
    std::function<real_t(real_t,real_t,real_t,real_t)> 
        generalFunction = [](real_t x, real_t y,real_t r,real_t t)
        {return 0.315*pow(abs(x),1.382)*pow(y,2.913)*pow(r,0.592)*pow(t,0.192);} ;

    real_t THRESHOLD = 2e-3;
    DREAM::FVM::fluxGridType fgTypes[4] = {DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, DREAM::FVM::FLUXGRIDTYPE_P1, DREAM::FVM::FLUXGRIDTYPE_P2, DREAM::FVM::FLUXGRIDTYPE_RADIAL};
    real_t maxError = 0;
    for(len_t fg  = 0; fg<4; fg++){
        len_t nr  = grid->GetNr() + (fg==3);
        len_t np1 = grid->GetNp1(0) + (fg==1);
        len_t np2 = grid->GetNp2(0) + (fg==2);
        for(len_t ir=0; ir<nr; ir++) 
            for(len_t i=0; i<np1; i++)
                for(len_t j=0; j<np2; j++){
                    real_t xi0, p;
                    if(fgTypes[fg] == DREAM::FVM::FLUXGRIDTYPE_P2){
                        xi0 = grid->GetMomentumGrid(0)->GetP2_f(j);
                        p = grid->GetMomentumGrid(0)->GetP1(i);
                    } else if(fgTypes[fg] == DREAM::FVM::FLUXGRIDTYPE_P1){
                        xi0 = grid->GetMomentumGrid(0)->GetP2(j);
                        p = grid->GetMomentumGrid(0)->GetP1_f(i);
                    } else { 
                        xi0 = grid->GetMomentumGrid(0)->GetP2(j);
                        p = grid->GetMomentumGrid(0)->GetP1(i);
                    }
                    real_t generalBounceAverage = grid->CalculateBounceAverage(ir,i,j,fgTypes[fg],generalFunction);
                    real_t generalBounceAverageAdaptive = grid_adaptive->CalculateBounceAverage(ir,i,j,fgTypes[fg],generalFunction);
                    real_t bounceAverageAtP = 0;//grid->GetRadialGrid()->CalculatePXiBounceAverageAtP(ir,xi0,fgTypes[fg],generalFunction);

                    real_t relativeError = 0;//getRelativeError(bounceAverageAtP, generalBounceAverage);
                    real_t relativeErrorAdaptive = 0;//getRelativeError(bounceAverageAtP, generalBounceAverageAdaptive);
                    real_t relativeErrorMethods  = getRelativeError(generalBounceAverageAdaptive, generalBounceAverage);

                    if(relativeError > maxError || relativeErrorAdaptive > maxError || relativeErrorMethods > maxError)
                        maxError = std::max(std::max(relativeError,relativeErrorAdaptive),relativeErrorMethods);
                        
                    if(!silentMode || relativeError >= THRESHOLD){
                        real_t r = grid->GetRadialGrid()->GetR(ir);
                        cout << "CompareBounceAverageMethods:" << endl;
                        cout << "----------------------------" << endl;
                        cout << "r: " << r  << ", xi0: " << xi0 << ", p: " << p << endl;
                        cout << "fluxGridType: " << fg << endl;
                        cout << "CalculateBounceAverage: " << generalBounceAverage << endl;
                        cout << "evaluateBounceAverageAtP: " << bounceAverageAtP << endl;
                        cout << "Relative error: " << relativeError << endl;
                        cout << "Relative error methods: " << relativeErrorMethods << endl << endl;
                        
                    }

                }
    }
    return maxError < THRESHOLD;
    
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
