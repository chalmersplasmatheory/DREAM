/**
 * Implementation of a class which handles the calculation of Coulomb logarithms.
 */
/**
 * The calculations are based on equations (2.7-10) from
 * L Hesslow et al., Generalized collision operator for fast electrons
 * interacting with partially ionized impurities, J Plasma Phys 84 (2018).
 */
#include "DREAM/Equations/CoulombLogarithm.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/FVMException.hpp"
#include <string>

using namespace DREAM;

/**
 * Construtor.
 */
CoulombLogarithm::CoulombLogarithm(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                enum OptionConstants::momentumgrid_type mgtype,  struct collqty_settings *cqset,
                len_t lnLambdaType):CollisionQuantity(g,u,ih,mgtype,cqset){
    // TODO: implement this in a better way, for now the Coulomb logarithm constructor
    // is called with an integer that determines whether it is an electron-electron lnLambda
    // or an electron-ion lnLambda.
    if(lnLambdaType==1)
        isLnEE = true;
    else if(lnLambdaType==2)
        isLnEI = true;
    else
        throw FVM::FVMException("Invalid lnLambdaType (supports only 1: electron-electron and 2: electron-ion)");
}

/**
 * Destructor.
 */
CoulombLogarithm::~CoulombLogarithm(){
    DeallocateCollisionQuantities();
    DeallocatePartialQuantities();
}

/**
 * Rebuilds the thermal Coulomb logarithm lnLambda_T and the relativistic
 * logarithm lnLambda_c (corresponding roughly to the energy-dependent one
 * evaluated at p=mc)
 */
void CoulombLogarithm::RebuildPlasmaDependentTerms(){
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    //real_t *n_cold = unknowns->GetUnknownData(id_ncold);
    real_t *n_free = new real_t[nr];
    n_free = ionHandler->evaluateFreeElectronDensityFromQuasiNeutrality(n_free);

    real_t *n_cold = unknowns->GetUnknownData(id_ncold);

    for(len_t ir=0; ir<nr; ir++){
        // should use n_free rather than n_cold; keeping it this way for now for benchmarking purposes
//        lnLambda_T[ir] = 14.9 - 0.5*log(n_free[ir]/1e20)  + log(T_cold[ir]/1e3);
        lnLambda_T[ir] = 14.9 - 0.5*log(n_cold[ir]/1e20)  + log(T_cold[ir]/1e3);
        lnLambda_c[ir] = 14.6 + 0.5*log( T_cold[ir]/(n_cold[ir]/1e20) );
        // should probably use the lower definition in the end, commented out for benchmarking 
        //lnLambda_c[ir] = lnLambda_T[ir] - 0.5*log(T_cold[ir]/Constants::mc2inEV);
    }
    delete [] n_free;
}

/**
 * Evaluates the Coulomb logarithm at radial index ir and momentum p.
 */
real_t CoulombLogarithm::evaluateAtP(len_t ir, real_t p){
    if(collQtySettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_CONSTANT){
        return lnLambda_c[ir];
    }
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    real_t gamma = sqrt(1+p*p);
    real_t eFactor;
    real_t pTeOverC = sqrt(2*T_cold[ir]/Constants::mc2inEV);
    if(isLnEE)
        eFactor = sqrt(2*(gamma-1))/pTeOverC;
    else if(isLnEI)
        eFactor = 2*p / pTeOverC;
    return lnLambda_T[ir] + log( 1 + pow(eFactor,kInterpolate) )/kInterpolate;
}


/**
 * Calculates and stores lnLambda. Calculates it in different ways depending on if the coulomb logarithm is the
 * electron-electron logarithm (isLnEE==true) or the electron-ion logarithm (isLnEI==true), or if we use 
 * the constant-logarithm approximation (lnL = lnL_c) or the full energy-dependent one.
 */
void CoulombLogarithm::AssembleQuantity(real_t **&collisionQuantity, len_t nr, len_t np1, len_t np2, len_t fluxGridType){
    const real_t *p;
    if(collQtySettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_CONSTANT){
        AssembleConstantLnLambda(collisionQuantity,nr,np1,np2);        
    } else if(isPXiGrid){
        // Optimized calculation for when a P-Xi grid is employed
        
        if(fluxGridType == 2)
            p = mg->GetP1_f();
        else
            p = mg->GetP1();        
        AssembleWithPXiGrid(collisionQuantity,p,nr,np1,np2);
    } else {
        if(fluxGridType == 2)
            p = mg->GetP_f1();
        else if(fluxGridType == 3)
            p = mg->GetP_f2();
        else
            p = mg->GetP();
        AssembleWithGeneralGrid(collisionQuantity, mg->GetP(),nr,np1,np2);
    }
}

/** 
 * Stores lnLambda when the constant lnLambda setting is used 
 * (then taking the relativistic value). 
 */
void CoulombLogarithm::AssembleConstantLnLambda(real_t **&lnLambda, len_t nr, len_t np1, len_t np2){
    for(len_t ir=0; ir<nr; ir++)
        for(len_t i=0; i<np1; i++)
            for(len_t j=0; j<np2; j++)
                lnLambda[ir][np1*j+i] = lnLambda_c[ir];
}

/**
 * Calculates and stores the energy-dependent lnLambda with a more efficient method when a P-Xi grid is utilized.
 */
void CoulombLogarithm::AssembleWithPXiGrid(real_t **&lnLambda,const real_t *pVec, len_t nr, len_t np1, len_t np2){
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    real_t p,gamma, pTeOverC, eFactor, lnL;
    for(len_t ir=0; ir<nr; ir++){
        pTeOverC = sqrt(2*T_cold[ir]/Constants::mc2inEV);
        for(len_t i = 0; i<np1; i++){
            p = pVec[i];
            gamma = sqrt(1+p*p);
            if(isLnEE)
                eFactor = sqrt(2*(gamma-1))/pTeOverC;
            else if(isLnEI)
                eFactor = 2*p / pTeOverC;

            lnL = lnLambda_T[ir] + log( 1 + pow(eFactor,kInterpolate) )/kInterpolate;
            for(len_t j=0; j<np2; j++){
                lnLambda[ir][np1*j+i] = lnL;
            }
        }
    }
}

/**
 * Calculates and stores the energy-dependent lnLambda when a non-pxi grid is used.
 */
void CoulombLogarithm::AssembleWithGeneralGrid(real_t **&lnLambda,const real_t *pVec, len_t nr, len_t np1, len_t np2){
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    real_t p,gamma, pTeOverC, eFactor;
    for(len_t ir=0; ir<nr; ir++){
        pTeOverC = sqrt(2*T_cold[ir]/Constants::mc2inEV);
        for(len_t i = 0; i<np1; i++){
            for(len_t j=0; j<np2; j++){
                p = pVec[np1*j+i];
                gamma = sqrt(1+p*p);
                if(isLnEE)
                    eFactor = sqrt(2*(gamma-1))/pTeOverC;
                else if(isLnEI)
                    eFactor = 2*p / pTeOverC;
                lnLambda[ir][np1*j+i] = lnLambda_T[ir] + log( 1 + pow(eFactor,kInterpolate) )/kInterpolate;
            }
        }
    }

}


void CoulombLogarithm::AllocatePartialQuantities(){
    DeallocatePartialQuantities();
    lnLambda_c = new real_t[nr];
    lnLambda_T = new real_t[nr];   
}

void CoulombLogarithm::DeallocatePartialQuantities(){
    if(lnLambda_c != nullptr)
        delete [] lnLambda_c;

    if(lnLambda_T != nullptr)
        delete [] lnLambda_T;
}
