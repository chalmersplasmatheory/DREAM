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
                CollisionQuantity::LnLambdaType lnLambdaType):CollisionQuantity(g,u,ih,mgtype,cqset){
    if(lnLambdaType==CollisionQuantity::LNLAMBDATYPE_EE)
        isLnEE = true;
    else if(lnLambdaType==CollisionQuantity::LNLAMBDATYPE_EI)
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
    for(len_t ir=0; ir<nr; ir++){
        lnLambda_T[ir] = evaluateLnLambdaT(ir);
        lnLambda_c[ir] = evaluateLnLambdaC(ir);
        lnLambda_ii[ir] = evaluateLnLambdaII(ir);
    }
}

/**
 * Rebuild method to be used in place of Rebuild() when only concerned about
 * using CoulombLogarithm::evaluateAtP(), such as in RunawayFluid.
 */
void CoulombLogarithm::RebuildRadialTerms(){
    if(gridRebuilt){
        nr  = rGrid->GetNr();
        AllocatePartialQuantities();
    }
    for(len_t ir=0; ir<nr; ir++){
        lnLambda_T[ir] = evaluateLnLambdaT(ir);
        lnLambda_c[ir] = evaluateLnLambdaC(ir);
        lnLambda_ii[ir] = evaluateLnLambdaII(ir);
    }
}


/**
 * Evaluates the relativistic lnLambda at radial grid index ir.
 */
real_t CoulombLogarithm::evaluateLnLambdaC(len_t ir){
    real_t T_cold = unknowns->GetUnknownData(id_Tcold)[ir];
    //real_t n_cold = unknowns->GetUnknownData(id_ncold)[ir];
    real_t n_free = ionHandler->GetFreeElectronDensityFromQuasiNeutrality(ir);
    if(n_free==0)
        return 0;
    return 14.6 + 0.5*log( T_cold/(n_free/1e20) );

    // eventually, to be more accurate, we may want to define it as
    // lnLambda_c = lnLambda_T - 0.5*log(T_cold/Constants::mc2inEV);
}


/**
 * Evaluates the thermal lnLambda at radial grid index ir.
 */
real_t CoulombLogarithm::evaluateLnLambdaT(len_t ir){
    real_t T_cold = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t n_free = ionHandler->GetFreeElectronDensityFromQuasiNeutrality(ir);
    if(n_free == 0)
        return 0;
    return 14.9 - 0.5*log(n_free/1e20) + log(T_cold/1e3);
}

/**
 * Evaluates the ion-ion lnLambda at radial grid index ir.
 */
real_t CoulombLogarithm::evaluateLnLambdaII(len_t ir){
    real_t T_cold = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t n_free = ionHandler->GetFreeElectronDensityFromQuasiNeutrality(ir);
    if(n_free == 0)
        return 0;
    return 17.3 - 0.5*log(n_free/1e20) + 1.5*log(T_cold/1e3);
}

/**
 * Evaluates the Coulomb logarithm at radial index ir and momentum p
 * using the provided collqty_settings.
 */
real_t CoulombLogarithm::evaluateAtP(len_t ir, real_t p,collqty_settings *inSettings){
    if(inSettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_CONSTANT)
        return  lnLambda_c[ir];
    else if(inSettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_THERMAL)
        return  lnLambda_T[ir];
    else if(inSettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_ION_ION)
        return  lnLambda_ii[ir];
    
    
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
 * Evaluates the Jacobian with respect to unknown derivId of the  
 * Coulomb logarithm at radial grid point ir and momentum p.
 */
real_t CoulombLogarithm::evaluatePartialAtP(len_t ir, real_t p, len_t derivId, len_t n,struct collqty_settings *inSettings){     
    if( (derivId != id_ni) && (derivId != id_Tcold) )
        return 0;
    
    if(derivId == id_ni){
        real_t n_free = ionHandler->GetFreeElectronDensityFromQuasiNeutrality(ir);
        if(n_free==0)
            return 0;
        len_t iz, Z0;
        ionHandler->GetIonIndices(n, iz, Z0);
        return -0.5 * Z0 / n_free; 
    }

    real_t T_cold = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t dlnL = 1/T_cold; // contribution from lnLambda_c
    if(inSettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_CONSTANT)
        return 0.5*dlnL;
    else if(inSettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_THERMAL)
        return dlnL;
    else if(inSettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_ION_ION)
        return 1.5*dlnL;
    

    // add the energy-dependent-lnLambda part    
    real_t eFactor;
    real_t pTeOverC = sqrt(2*T_cold/Constants::mc2inEV);
    if(isLnEE) {
        real_t gamma = sqrt(1+p*p);
        eFactor = sqrt(2*(gamma-1))/pTeOverC;
    } else if(isLnEI)
        eFactor = 2*p / pTeOverC;

    // derivative of eFactor wrt Tcold
    real_t deFac = - eFactor / (2*T_cold); 
    real_t powFac = pow(eFactor,kInterpolate-1);
    dlnL += powFac * deFac/(1 + powFac*eFactor); 
    
    return dlnL;
}


/**
 * Calculates and stores lnLambda. Calculates it in different ways depending on if the coulomb logarithm is the
 * electron-electron logarithm (isLnEE==true) or the electron-ion logarithm (isLnEI==true), or if we use 
 * the constant-logarithm approximation (lnL = lnL_c) or the full energy-dependent one.
 */
void CoulombLogarithm::AssembleQuantity(real_t **&collisionQuantity, len_t nr, len_t np1, len_t np2, FVM::fluxGridType fluxGridType){
    const real_t *p;
    if( (collQtySettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_CONSTANT)
     || (collQtySettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_THERMAL) 
     || (collQtySettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_ION_ION) ){
        AssembleConstantLnLambda(collisionQuantity,nr,np1,np2);        
    } else if(isPXiGrid){
        // Optimized calculation for when a P-Xi grid is employed
        if(fluxGridType == FVM::FLUXGRIDTYPE_P1)
            p = mg->GetP1_f();
        else
            p = mg->GetP1();        
        AssembleWithPXiGrid(collisionQuantity,p,nr,np1,np2);
    } else {
        if(fluxGridType == FVM::FLUXGRIDTYPE_P1)
            p = mg->GetP_f1();
        else if(fluxGridType == FVM::FLUXGRIDTYPE_P2)
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
    for(len_t ir=0; ir<nr; ir++){
        real_t lnL;
        if(collQtySettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_CONSTANT)
            lnL =  lnLambda_c[ir];
        else if (collQtySettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_THERMAL)
            lnL = lnLambda_T[ir];
        else if (collQtySettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_ION_ION)
            lnL = lnLambda_ii[ir];
        for(len_t i=0; i<np1*np2; i++)
            lnLambda[ir][i] = lnL;
    }
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
            for(len_t j=0; j<np2; j++)
                lnLambda[ir][np1*j+i] = lnL;
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
        for(len_t i = 0; i<np1; i++)
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


/**
 * Allocate quantities
 */
void CoulombLogarithm::AllocatePartialQuantities(){
    DeallocatePartialQuantities();
    lnLambda_c = new real_t[nr];
    lnLambda_T = new real_t[nr];  
    lnLambda_ii = new real_t[nr]; 
}


/**
 * Deallocate quantities
 */
void CoulombLogarithm::DeallocatePartialQuantities(){
    if(lnLambda_c != nullptr)
        delete [] lnLambda_c;
    if(lnLambda_T != nullptr)
        delete [] lnLambda_T;
    if(lnLambda_ii != nullptr)
        delete [] lnLambda_ii;
}
