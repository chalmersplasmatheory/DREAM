#include "DREAM/Equations/CoulombLogarithm.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/FVMException.hpp"
#include <string>

using namespace DREAM;

CoulombLogarithm::CoulombLogarithm(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantityHandler::collqtyhand_settings *cqset,
                len_t lnLambdaType):CollisionQuantity(g,u,ih,mgtype,cqset){
    if(lnLambdaType==1)
        isLnEE = true;
    else if(lnLambdaType==2)
        isLnEI = true;
    else
        throw FVM::FVMException("Invalid lnLambdaType (supports only 1: electron-electron and 2: electron-ion)");
}

CoulombLogarithm::~CoulombLogarithm(){
    DeallocatePartialQuantities();
}

void CoulombLogarithm::RebuildPlasmaDependentTerms(){
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    //real_t *n_cold = unknowns->GetUnknownData(id_ncold);
    real_t *n_free = new real_t[nr];
    n_free = ionHandler->evaluateFreeElectronDensityFromQuasiNeutrality(n_free);

    for(len_t ir=0; ir<nr; ir++){
        lnLambda_T[ir] = 14.9 - 0.5*log(n_free[ir]/1e20)  + log(T_cold[ir]/1e3);
        lnLambda_c[ir] = lnLambda_T[ir] - 0.5*log(T_cold[ir]/Constants::mc2inEV);
    }
    delete [] n_free;
}
/*
        lnLTe[ir] = 14.9 + 0.5*log( (T_cold[ir]/1e3)*(T_cold[ir]/1e3)/(n_cold[ir]/1e20) );
        lnLc[ir] = lnLTe[ir] - 0.5*log(T_cold[ir]/Constants::mc2inEV);
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
 * the constant-logarithm approximation (lnL = lnL_c). 
 */
void CoulombLogarithm::AssembleQuantity(){

    if(collQtySettings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_CONSTANT){
        AssembleConstantLnLambda();        
    } else if(isPXiGrid)
        // Optimized calculation for when a P-Xi grid is employed
        AssembleWithPXiGrid();
    else 
        AssembleWithGeneralGrid();
}


void CoulombLogarithm::AssembleConstantLnLambda(){
    if(!buildOnlyF1F2){
        AssembleConstantLnLambda(collisionQuantity,nr,np1,np2);
        AssembleConstantLnLambda(collisionQuantity_fr,nr+1,np1,np2);
    }
    AssembleConstantLnLambda(collisionQuantity_f1,nr,np1+1,np2);
    AssembleConstantLnLambda(collisionQuantity_f2,nr,np1,np2+1);
        
}
void CoulombLogarithm::AssembleConstantLnLambda(real_t **&lnLambda, len_t nr, len_t np1, len_t np2){
    for(len_t ir=0; ir<nr; ir++)
        for(len_t i=0; i<np1; i++)
            for(len_t j=0; j<np2; j++)
                lnLambda[ir][np1*j+i] = lnLambda_c[ir];
}


void CoulombLogarithm::AssembleWithPXiGrid(){
    if(!buildOnlyF1F2){
        AssembleWithPXiGrid(collisionQuantity, mg->GetP1(),nr,np1,np2);
        AssembleWithPXiGrid(collisionQuantity_fr,mg->GetP1(),nr+1,np1,np2);
    }
    AssembleWithPXiGrid(collisionQuantity_f1,mg->GetP1_f(),nr,np1+1,np2);
    AssembleWithPXiGrid(collisionQuantity_f2,mg->GetP1(),nr,np1,np2+1);

}

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


void CoulombLogarithm::AssembleWithGeneralGrid(){
    if(!buildOnlyF1F2){
        AssembleWithGeneralGrid(collisionQuantity, mg->GetP(),nr,np1,np2);
        AssembleWithGeneralGrid(collisionQuantity_fr,mg->GetP(),nr+1,np1,np2);
    }
    AssembleWithGeneralGrid(collisionQuantity_f1,mg->GetP_f1(),nr,np1+1,np2);
    AssembleWithGeneralGrid(collisionQuantity_f2,mg->GetP_f2(),nr,np1,np2+1);

}


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
