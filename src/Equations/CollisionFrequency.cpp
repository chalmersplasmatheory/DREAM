/**
 * Implementation of the a which handles the construction of 
 * collision frequencies as well as their partial derivatives 
 * with respect to the unknown quantities.
 * CollisionFrequencies should be rebuilt after their CoulombLogarithms
 * have been rebuilt.
 */

#include "DREAM/Equations/CollisionFrequency.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/FVMException.hpp"
#include <string>
#include "gsl/gsl_sf_bessel.h"

using namespace DREAM;

/**
 * Constructor.
 */
CollisionFrequency::CollisionFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                CoulombLogarithm *lnLee, CoulombLogarithm *lnLei,
                enum OptionConstants::momentumgrid_type mgtype,  struct collqty_settings *cqset)
                : CollisionQuantity(g,u,ih,mgtype,cqset) {
    lnLambdaEE = lnLee;
    lnLambdaEI = lnLei;
}


/**
 * Destructor.
 */
CollisionFrequency::~CollisionFrequency(){
    DeallocatePartialQuantities();
}


/**
 * Evaluates the collision frequency at radial grid point ir and momentum p,
 * neglecting any contribution from the nonlinear collision operator, using
 * the input collqty_settings object.
 */
real_t CollisionFrequency::evaluateAtP(len_t ir, real_t p,collqty_settings *inSettings){ 
    bool isPartiallyScreened = (inSettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED);
    bool isNonScreened = (inSettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED);
    bool isBrems = (inSettings->bremsstrahlung_mode != OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT);
    real_t *ncold = unknowns->GetUnknownData(id_ncold);
    real_t ntarget = ncold[ir];
    if (isNonScreened)
        ntarget += ionHandler->GetBoundElectronDensity(ir);

    real_t preFact = evaluatePreFactorAtP(p,inSettings->collfreq_mode); 
    real_t lnLee = lnLambdaEE->evaluateAtP(ir,p,inSettings);
    real_t lnLei = lnLambdaEI->evaluateAtP(ir,p,inSettings);
    
    // Add electron contribution to collision frequency
    real_t collFreq = lnLee * evaluateElectronTermAtP(ir,p,inSettings->collfreq_mode) * ntarget;

    len_t ind;
    // Add ion contribution; SlowingDownFrequency doesn't have one and will skip this step
    if(hasIonTerm){
        for(len_t iz = 0; iz<nZ; iz++)
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                ind = ionIndex[iz][Z0];
                len_t Zfact = Z0*Z0;
                if(isNonScreened)
                    Zfact = Zs[iz]*Zs[iz];
                collFreq += lnLei * Zfact * evaluateIonTermAtP(iz,Z0,p) * ionDensities[ir][ind];
            }
    }
    // Add screening contribution
    if(isPartiallyScreened)
        for(len_t iz = 0; iz<nZ; iz++)
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                ind = ionIndex[iz][Z0];
                collFreq +=  evaluateScreenedTermAtP(iz,Z0,p,inSettings->collfreq_mode) * ionDensities[ir][ind];
            }
    collFreq *= preFact;

    // Add Bremsstrahlung contribution
    if(isBrems)
        for(len_t iz = 0; iz<nZ; iz++)
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                ind = ionIndex[iz][Z0];
                collFreq +=  evaluateBremsstrahlungTermAtP(iz,Z0,p,inSettings->bremsstrahlung_mode,inSettings->collfreq_type) * ionDensities[ir][ind];
            }

    return collFreq;
}


/**
 *  Calculates and stores partial terms which depend on unknown quantities (density and temperature).
 */
void CollisionFrequency::RebuildPlasmaDependentTerms(){
    len_t indZ;
    for(len_t iz = 0; iz<nZ; iz++)
        for(len_t Z0=0; Z0<=Zs[iz]; Z0++)
            for(len_t ir=0; ir<nr; ir++){
                indZ = ionIndex[iz][Z0];            
                ionDensities[ir][indZ] = ionHandler->GetIonDensity(ir,iz,Z0);
            }
    
    if(collQtySettings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL)
        InitializeGSLWorkspace();
    if (!buildOnlyF1F2){
        setElectronTerm(nColdTerm,mg->GetP(),nr,np1,np2_store);
        setElectronTerm(nColdTerm_fr,mg->GetP(),nr/*+1*/,np1,np2_store);
    }
    setElectronTerm(nColdTerm_f1,mg->GetP_f1(),nr,np1+1,np2_store);
    setElectronTerm(nColdTerm_f2,mg->GetP_f2(),nr,np1,np2_store+1);
}


/**
 * Calculates and stores all quantities that do not depend on momentum. Is to be used instead of Rebuild when
 * you only want to use the CollisionFrequency::evaluateAtP function (such as in RunawayFluid).
 */
void CollisionFrequency::RebuildRadialTerms(){
    if(gridRebuilt){
        DeallocateRadialQuantities();
        nr  = rGrid->GetNr();
        nZ  = ionHandler->GetNZ();
        nzs = ionHandler->GetNzs();
        AllocateRadialQuantities();
    }
    len_t indZ;
    const len_t *ZAtomicCharge = ionHandler->GetZs();
    for(len_t iz = 0; iz<nZ; iz++){
        Zs[iz] = ZAtomicCharge[iz];
        for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
            indZ = ionHandler->GetIndex(iz,Z0);
            ionIndex[iz][Z0] = indZ; 
            atomicParameter[indZ] = GetAtomicParameter(iz,Z0);
            for(len_t ir=0; ir<nr; ir++)
                ionDensities[ir][indZ] = ionHandler->GetIonDensity(ir,iz,Z0);
        }
    }
    if(collQtySettings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL)
        InitializeGSLWorkspace();
}


/**
 * Allocates quantities involved in RebuildRadialTerms()
 */
void CollisionFrequency::AllocateRadialQuantities(){    
//    DeallocateRadialQuantities();
    InitializeGSLWorkspace();
    Zs = new real_t[nZ];
    ionIndex = new real_t*[nZ];
    ionDensities = new real_t*[nr];
    atomicParameter = new real_t[nzs];

    for(len_t iz=0;iz<nZ;iz++)
        ionIndex[iz] = new real_t[ionHandler->GetZ(iz)+1];
    for(len_t ir=0; ir<nr;ir++)
        ionDensities[ir] = new real_t[nzs];

}


/**
 * Deallocates quantities involved in RebuildRadialTerms()
 */
void CollisionFrequency::DeallocateRadialQuantities(){
    if(Zs!=nullptr){
        delete [] Zs;
        delete [] atomicParameter;

        for(len_t iz=0;iz<nZ;iz++)
            delete [] ionIndex[iz];
        for(len_t ir=0; ir<nr;ir++)
            delete [] ionDensities[ir];
        delete [] ionIndex;
        delete [] ionDensities;
    }
}


/**
 * Calculates and stores partial contributions that only depend on the grid. If using P-Xi grid, 
 * only store momentum dependent quantities on a size np1 array. 
 */
void CollisionFrequency::RebuildConstantTerms(){
    const len_t *ZAtomicCharge = ionHandler->GetZs();
    len_t indZ;
    for(len_t iz = 0; iz<nZ; iz++){
        Zs[iz] = ZAtomicCharge[iz];
        for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
            indZ = ionHandler->GetIndex(iz,Z0);
            ionIndex[iz][Z0] = indZ; 
            atomicParameter[indZ] = GetAtomicParameter(iz,Z0);
        }
    }

    if (!buildOnlyF1F2){
        setPreFactor(preFactor,mg->GetP(),np1,np2_store);
        setPreFactor(preFactor_fr,mg->GetP(),np1,np2_store);
        setIonTerm(ionTerm,mg->GetP(),np1,np2_store);
        setIonTerm(ionTerm_fr,mg->GetP(),np1,np2_store);
        if(isBrems){
            setBremsTerm(bremsTerm,mg->GetP(),np1,np2_store);
            setBremsTerm(bremsTerm_fr,mg->GetP(),np1,np2_store);
        }
        if(isPartiallyScreened){
            setScreenedTerm(screenedTerm,mg->GetP(),np1,np2_store);
            setScreenedTerm(screenedTerm_fr,mg->GetP(),np1,np2_store);
        }
    }
    setPreFactor(preFactor_f1,mg->GetP_f1(),np1+1,np2_store);
    setPreFactor(preFactor_f2,mg->GetP_f2(),np1,np2_store+1);
    setIonTerm(ionTerm_f1,mg->GetP_f1(),np1+1,np2_store);
    setIonTerm(ionTerm_f2,mg->GetP_f2(),np1,np2_store+1);
    if(isBrems){
        setBremsTerm(bremsTerm_f1,mg->GetP_f1(),np1+1,np2_store);
        setBremsTerm(bremsTerm_f2,mg->GetP_f2(),np1,np2_store+1);
    }
    if(isPartiallyScreened){
        setScreenedTerm(screenedTerm_f1,mg->GetP_f1(),np1+1,np2_store);
        setScreenedTerm(screenedTerm_f2,mg->GetP_f2(),np1,np2_store+1);
    }
    if(isNonlinear)
        calculateIsotropicNonlinearOperatorMatrix();
}


/**
 * Calculates and stores the partial contributions (i.e. kind of partial derivates of 
 * the collision frequencies)
 */
void CollisionFrequency::SetPartialContributions(FVM::fluxGridType fluxGridType){
    if(fluxGridType==FVM::FLUXGRIDTYPE_DISTRIBUTION){
        SetNColdPartialContribution(nColdTerm,preFactor,lnLambdaEE->GetValue(),nr,np1,np2,nColdPartialContribution);
        SetNiPartialContribution(nColdTerm,ionTerm, screenedTerm,bremsTerm,preFactor,lnLambdaEE->GetValue(),lnLambdaEI->GetValue(),nr,np1,np2,ionPartialContribution, ionLnLambdaPartialContribution);
        SetTColdPartialContribution(nColdTerm,ionTerm, preFactor,lnLambdaEE->GetValue(),mg->GetP(), nr,np1,np2,TColdPartialContribution);
    } else if(fluxGridType==FVM::FLUXGRIDTYPE_RADIAL){
        SetNColdPartialContribution(nColdTerm_fr,preFactor_fr,lnLambdaEE->GetValue_fr(),nr /*+1*/,np1,np2,nColdPartialContribution_fr);
        SetNiPartialContribution(nColdTerm_fr,ionTerm_fr,screenedTerm_fr,bremsTerm_fr, preFactor_fr,lnLambdaEE->GetValue_fr(),lnLambdaEI->GetValue_fr(),nr/*+1*/,np1,np2,ionPartialContribution_fr, ionLnLambdaPartialContribution_fr);
        SetTColdPartialContribution(nColdTerm_fr,ionTerm_fr,preFactor_fr,lnLambdaEE->GetValue_fr(),mg->GetP(), nr/*+1*/,np1,np2,TColdPartialContribution_fr);
    } else if(fluxGridType==FVM::FLUXGRIDTYPE_P1){
        SetNColdPartialContribution(nColdTerm_f1,preFactor_f1,lnLambdaEE->GetValue_f1(),nr,np1+1,np2,nColdPartialContribution_f1);
        SetNiPartialContribution(nColdTerm_f1,ionTerm_f1,screenedTerm_f1,bremsTerm_f1, preFactor_f1,lnLambdaEE->GetValue_f1(),lnLambdaEI->GetValue_f1(),nr,np1+1,np2,ionPartialContribution_f1, ionLnLambdaPartialContribution_f1);
        SetTColdPartialContribution(nColdTerm_f1,ionTerm_f1,preFactor_f1,lnLambdaEE->GetValue_f1(),mg->GetP_f1(), nr,np1+1,np2,TColdPartialContribution_f1);
    } else if(fluxGridType==FVM::FLUXGRIDTYPE_P2){
        SetNColdPartialContribution(nColdTerm_f2,preFactor_f2,lnLambdaEE->GetValue_f2(),nr,np1,np2+1,nColdPartialContribution_f2);
        SetNiPartialContribution(nColdTerm_f2,ionTerm_f2,screenedTerm_f2,bremsTerm_f2, preFactor_f2,lnLambdaEE->GetValue_f2(),lnLambdaEI->GetValue_f2(),nr,np1,np2+1,ionPartialContribution_f2, ionLnLambdaPartialContribution_f2);
        SetTColdPartialContribution(nColdTerm_f2,ionTerm_f2,preFactor_f2,lnLambdaEE->GetValue_f2(),mg->GetP_f2(), nr,np1,np2+1,TColdPartialContribution_f2);
    }
    if(isNonlinear && (fluxGridType == FVM::FLUXGRIDTYPE_P1) )
        SetNonlinearPartialContribution(lnLambdaEE,fHotPartialContribution_f1);
}


/**
 * Calculates and stores the partial contributions to the collision frequency, 
 * and puts them together to get the full thing.
 */
void CollisionFrequency::AssembleQuantity(real_t **&collisionQuantity,  len_t nr, len_t np1, len_t np2, enum FVM::fluxGridType fluxGridType){
    real_t collQty;
    real_t *ncold = unknowns->GetUnknownData(id_ncold);

    SetPartialContributions(fluxGridType);

    const real_t *nColdContribution = GetNColdPartialContribution(fluxGridType);
    real_t *ionLnLContrib;
    const real_t *ionContribution = GetNiPartialContribution(fluxGridType, &ionLnLContrib);

    for(len_t ir=0; ir<nr; ir++)
        for(len_t pind=0; pind<np1*np2; pind++){
            // the collision frequencies are linear in ncold
            collQty = ncold[ir]*nColdContribution[np1*np2*ir + pind];
            for(len_t indZ = 0; indZ<nzs; indZ++){
                len_t ind = (indZ*nr + ir)*np1*np2 + pind;
                // when subtracting the lnLambda terms, the collision frequencies are linear in ion densities
                collQty += ionDensities[ir][indZ]*(ionContribution[ind] - ionLnLContrib[ind]);
            }
            collisionQuantity[ir][pind] = collQty; 
        }
}


/**
 * Calculation of the partial contribution to the collision frequency from the unknown quantity
 * with ID id_unknown. Returns the partial derivative of the term with respect to that quantity 
 * (ignoring variations with lnLambda). See AssembleQuantity(...) or addNonlinearContribution() 
 * for how it is used.
 */
const real_t* CollisionFrequency::GetUnknownPartialContribution(len_t id_unknown, FVM::fluxGridType fluxGridType) const{
    if(id_unknown == id_ncold)
        return GetNColdPartialContribution(fluxGridType);
    else if(id_unknown == id_ni)
        return GetNiPartialContribution(fluxGridType);
    else if(id_unknown == id_Tcold)
        return GetTColdPartialContribution(fluxGridType);
    else if(id_unknown == unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT)){
        if(!( (fluxGridType==FVM::FLUXGRIDTYPE_P1)&&(np2==1)&&(isPXiGrid) ) )
            throw FVM::FVMException("Nonlinear contribution to collision frequencies is only implemented for hot-tails, with p-xi grid and np2=1 and evaluated on the p flux grid.");
        return GetNonlinearPartialContribution(fluxGridType);
    } else {
        return nullptr;
//        throw FVM::FVMException("Invalid id_unknown: %s does not contribute to the collision frequencies",unknowns->GetUnknown(id_unknown)->GetName());
    }
}

const real_t* CollisionFrequency::GetNColdPartialContribution(FVM::fluxGridType fluxGridType) const{
    if(fluxGridType==FVM::FLUXGRIDTYPE_DISTRIBUTION)
        return nColdPartialContribution;
    else if (fluxGridType==FVM::FLUXGRIDTYPE_RADIAL)
        return nColdPartialContribution_fr;
    else if (fluxGridType==FVM::FLUXGRIDTYPE_P1)
        return nColdPartialContribution_f1;
    else if (fluxGridType==FVM::FLUXGRIDTYPE_P2)
        return nColdPartialContribution_f2;
    else {
        throw FVM::FVMException("Invalid fluxGridType");
        return nullptr;
    }
}

const real_t* CollisionFrequency::GetNiPartialContribution(FVM::fluxGridType fluxGridType, real_t **lnLambdaContrib) const{
    bool setLL = (lnLambdaContrib != nullptr);
    if(fluxGridType==FVM::FLUXGRIDTYPE_DISTRIBUTION){
        if(setLL) *lnLambdaContrib = ionLnLambdaPartialContribution;
        return ionPartialContribution;
    } else if (fluxGridType==FVM::FLUXGRIDTYPE_RADIAL){
        if(setLL) *lnLambdaContrib = ionLnLambdaPartialContribution_fr;
        return ionPartialContribution_fr;
    } else if (fluxGridType==FVM::FLUXGRIDTYPE_P1){
        if(setLL) *lnLambdaContrib = ionLnLambdaPartialContribution_f1;
        return ionPartialContribution_f1;
    } else if (fluxGridType==FVM::FLUXGRIDTYPE_P2){
        if(setLL) *lnLambdaContrib = ionLnLambdaPartialContribution_f2;
        return ionPartialContribution_f2;
    } else {
        throw FVM::FVMException("Invalid fluxGridType");
        return nullptr;
    }
}

const real_t* CollisionFrequency::GetTColdPartialContribution(FVM::fluxGridType fluxGridType) const{
    if(fluxGridType==FVM::FLUXGRIDTYPE_DISTRIBUTION)
        return TColdPartialContribution;
    else if (fluxGridType==FVM::FLUXGRIDTYPE_RADIAL)
        return TColdPartialContribution_fr;
    else if (fluxGridType==FVM::FLUXGRIDTYPE_P1)
        return TColdPartialContribution_f1;
    else if (fluxGridType==FVM::FLUXGRIDTYPE_P2)
        return TColdPartialContribution_f2;
    else {
        throw FVM::FVMException("Invalid fluxGridType");
        return nullptr;
    }
}

const real_t* CollisionFrequency::GetNonlinearPartialContribution(FVM::fluxGridType fluxGridType) const{
    if(fluxGridType==FVM::FLUXGRIDTYPE_P1)
        return fHotPartialContribution_f1;
    else {
//        throw FVM::FVMException("Invalid fluxGridType. Nonlinear contribution only supported for p1 flux grid.");
        return nullptr;
    }
}


/** Adds the non-linear contribution to the collision frequency. For now, only supports 
 * hot-tail grids where np2=1 and using a pxi-grid, and only updates the p flux grid 
 * component.
 */
void CollisionFrequency::AddNonlinearContribution(){
    real_t *fHot = unknowns->GetUnknownData(OptionConstants::UQTY_F_HOT);
    const real_t* const fHotPartialContribution_f1 = GetNonlinearPartialContribution(FVM::FLUXGRIDTYPE_P1);

    for (len_t ir=0;ir<nr;ir++)
        for(len_t i=0; i<np1+1; i++)
            for(len_t ip=0; ip<np1; ip++)
                collisionQuantity_f1[ir][i] += fHotPartialContribution_f1[ip*(np1+1)*nr + ir*(np1+1) + i] * fHot[np1*ir+ip];
}


/**
 * Calculates and stores the momentum-dependent prefactor to the collision frequencies.
 */
void CollisionFrequency::setPreFactor(real_t *&preFactor, const real_t *pIn, len_t np1, len_t np2){
    for (len_t pind = 0; pind<np1*np2; pind++)
        preFactor[pind] = evaluatePreFactorAtP(pIn[pind],collQtySettings->collfreq_mode);
}


/**
 * Calculates and stores the ion contribution to the collision frequency.
 */
void CollisionFrequency::setIonTerm(real_t *&ionTerm, const real_t *pIn, len_t np1, len_t np2){
    if(!hasIonTerm)
        return;
    real_t p;
    len_t ind, pind;
    for(len_t i = 0; i<np1; i++)
        for (len_t j = 0; j<np2; j++){
            pind = np1*j+i;
            p = pIn[pind];
            for(len_t iz = 0; iz<nZ; iz++)
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                    ind = ionIndex[iz][Z0];
                    ionTerm[ind*np1*np2 + pind] = evaluateIonTermAtP(iz,Z0,p);
                }
        }
}


/**
 * Calculates and stores the partially-screened contribution to the collision frequency.
 */
void CollisionFrequency::setScreenedTerm(real_t *&screenedTerm, const real_t *pIn, len_t np1, len_t np2){
    len_t ind;
    if(isPXiGrid)
        for(len_t i = 0; i<np1; i++)
            for(len_t iz = 0; iz<nZ; iz++)
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                    ind = ionIndex[iz][Z0];
                    real_t screenedAtP = evaluateScreenedTermAtP(iz,Z0,pIn[i], collQtySettings->collfreq_mode);
                    for (len_t j = 0; j<np2; j++)
                        screenedTerm[ind*np1*np2 + np1*j + i] = screenedAtP;
                }
    else
        for (len_t pind = 0; pind<np1*np2; pind++)
            for(len_t iz = 0; iz<nZ; iz++)
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                    ind = ionIndex[iz][Z0];
                    screenedTerm[ind*np1*np2 + pind] = evaluateScreenedTermAtP(iz,Z0,pIn[pind], collQtySettings->collfreq_mode);
                }
}


/**
 * Calculates and stores the bremsstrahlung contribution to the collision frequency.
 */
void CollisionFrequency::setBremsTerm(real_t *&bremsTerm, const real_t *pIn, len_t np1, len_t np2){
    len_t ind;

    if(isPXiGrid){
        for(len_t i = 0; i<np1; i++)
            for(len_t iz = 0; iz<nZ; iz++)
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                    ind = ionIndex[iz][Z0];
                    real_t bremsAtP = evaluateBremsstrahlungTermAtP(iz, Z0, pIn[i], collQtySettings->bremsstrahlung_mode, collQtySettings->collfreq_type);
                    for (len_t j = 0; j<np2; j++)
                        bremsTerm[ind*np1*np2 + np1*j + i] = bremsAtP;
                }
    } else {
        for (len_t pind = 0; pind<np1*np2; pind++)
            for(len_t iz = 0; iz<nZ; iz++)
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                    ind = ionIndex[iz][Z0];
                    bremsTerm[ind*np1*np2 + pind] = evaluateBremsstrahlungTermAtP(iz, Z0, pIn[pind], collQtySettings->bremsstrahlung_mode, collQtySettings->collfreq_type);
                }
    }
}


/**
 * Calculates and stores the free-electron contribution to the collision frequency.
 */
void CollisionFrequency::setElectronTerm(real_t **&nColdTerm, const real_t *pIn, len_t nr, len_t np1, len_t np2){
    if(isPXiGrid)
        for(len_t i=0;i<np1;i++)
            for(len_t ir=0; ir<nr; ir++){
                real_t electronTerm = evaluateElectronTermAtP(ir,pIn[i],collQtySettings->collfreq_mode);
                for(len_t j=0;j<np2;j++)
                    nColdTerm[ir][np1*j+i] = electronTerm;
            }
    else
        for(len_t pind=0; pind<np1*np2;pind++)
            for(len_t ir=0; ir<nr; ir++)
                nColdTerm[ir][pind] = evaluateElectronTermAtP(ir,pIn[pind],collQtySettings->collfreq_mode);
            
}


// PSI FUNCTIONS FOR EVALUATION OF "FULL" COLLFREQ_MODE 
// Is the relativistic generalisation of the Chandrasekhar functions,
// appearing in all collision frequencies.
real_t CollisionFrequency::psi0Integrand(real_t s, void *params){
    real_t Theta = *(real_t *) params;
    real_t gs = sqrt(1+s*s);
    real_t gsMinusOne = s*s/(1+gs); // = gs - 1
    return exp(-gsMinusOne/Theta)/gs;
}

real_t CollisionFrequency::psi1Integrand(real_t s, void *params){
    real_t Theta = *(real_t *) params;
    real_t gs = sqrt(1+s*s);
    real_t gsMinusOne = s*s/(1+gs); // = gs - 1
    return exp(-gsMinusOne/Theta);
}
real_t CollisionFrequency::psi2Integrand(real_t s, void *params){
    real_t Theta = *(real_t *) params;
    real_t gs = sqrt(1+s*s);
    real_t gsMinusOne = s*s/(1+gs); // = gs - 1
    return exp(-gsMinusOne/Theta)*gs;
}

/**
 * Evaluates the Psi0 thermal collision frequency function.
 */
real_t CollisionFrequency::evaluatePsi0(len_t ir, real_t p) {
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    real_t Theta = T_cold[ir] / Constants::mc2inEV;

    gsl_function F;
    F.function = &(CollisionFrequency::psi0Integrand); 
    F.params = &Theta;
    real_t psi0int, error; 

    real_t epsabs = 0, epsrel = 5e-4, lim = gsl_ad_w->limit; 
    gsl_integration_qag(&F,0,p,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w,&psi0int,&error);
    return psi0int;
}

/**
 * Evaluates the Psi1 thermal collision frequency function.
 */
real_t CollisionFrequency::evaluatePsi1(len_t ir, real_t p) {
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    real_t Theta = T_cold[ir] / Constants::mc2inEV;

    gsl_function F;
    F.function = &(CollisionFrequency::psi1Integrand); 
    F.params = &Theta;
    real_t psi1int, error; 

    real_t epsabs = 0, epsrel = 5e-4, lim = gsl_ad_w->limit; 
    gsl_integration_qag(&F,0,p,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w,&psi1int,&error);
    return psi1int;
}

real_t CollisionFrequency::evaluatePsi2(len_t ir, real_t p) {
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    real_t Theta = T_cold[ir] / Constants::mc2inEV;

    gsl_function F;
    F.function = &(CollisionFrequency::psi2Integrand); 
    F.params = &Theta;
    real_t psi2int, error; 

    real_t epsabs = 0, epsrel = 5e-4, lim = gsl_ad_w->limit; 
    gsl_integration_qags(&F,0,p,epsabs,epsrel,lim,gsl_ad_w,&psi2int,&error);
    return psi2int;
}


/**
 * Evaluates e^x K_n(x), with K_n the exponentially decreasing modified bessel function.
 */
real_t CollisionFrequency::evaluateExp1OverThetaK(real_t Theta, real_t n) {
    return gsl_sf_bessel_Kn_scaled(n,1.0/Theta);
}


void CollisionFrequency::SetNiPartialContribution(real_t **nColdTerm, real_t *ionTerm, real_t *screenedTerm, real_t *bremsTerm, real_t *preFactor, real_t *const* lnLee,  real_t *const* lnLei, len_t nr, len_t np1, len_t np2, real_t *&partQty, real_t *&ionLnLContrib){
    if(partQty==nullptr){
        partQty = new real_t[nzs*np1*np2*nr];
        ionLnLContrib = new real_t[nzs*np1*np2*nr];
    }
    for(len_t it = 0; it<nzs*np1*np2*nr; it++){
        partQty[it] = 0;
        ionLnLContrib[it] = 0;
    }

    real_t **lnLEE_partialNi = new real_t*[nr];
    real_t **lnLEI_partialNi = new real_t*[nr];
    for(len_t ir=0; ir<nr; ir++){
        lnLEE_partialNi[ir] = new real_t[nzs];
        lnLEI_partialNi[ir] = new real_t[nzs];
        for(len_t iz=0; iz<nZ; iz++)
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                len_t indZ = ionIndex[iz][Z0]; 
                lnLEE_partialNi[ir][indZ] = lnLambdaEE->evaluatePartialAtP(ir,0,id_ni,indZ);
                lnLEI_partialNi[ir][indZ] = lnLambdaEI->evaluatePartialAtP(ir,0,id_ni,indZ);
            }
    }
    len_t pind, indZ;
    real_t partContrib;
    real_t electronTerm;

    len_t pindStore;

    if(isPXiGrid) {
        for(len_t ir = 0; ir<nr; ir++){
            real_t ntarget = unknowns->GetUnknownData(id_ncold)[ir];
            if (isNonScreened)
                ntarget += ionHandler->GetBoundElectronDensity(ir);
            for(len_t i = 0; i<np1; i++){
                electronTerm = ntarget*nColdTerm[ir][i]*preFactor[i];
                for(len_t indZ=0; indZ<nzs; indZ++){
                    real_t lnLContrib = electronTerm * lnLEE_partialNi[ir][indZ];
                    for(len_t j = 0; j<np2; j++){
                        pind = np1*j+i;
                        len_t ind = (indZ*nr+ir)*np1*np2 + pind;
                        ionLnLContrib[ind] += lnLContrib; 
                        partQty[ind] += lnLContrib;
                    }
                }
            }      
        }
    } else {
        for(len_t i = 0; i<np1; i++)
            for(len_t j = 0; j<np2; j++){
                pind = np1*j+i;                
                for(len_t ir = 0; ir<nr; ir++){
                    real_t ntarget = unknowns->GetUnknownData(id_ncold)[ir];
                    if (isNonScreened)
                        ntarget += ionHandler->GetBoundElectronDensity(ir);
                    electronTerm = ntarget*nColdTerm[ir][pind]*preFactor[pind];
                    for(len_t indZ=0; indZ<nzs; indZ++){
                        len_t ind = (indZ*nr+ir)*np1*np2 + pind;
                        real_t lnLContrib = electronTerm * lnLEE_partialNi[ir][indZ];
                        ionLnLContrib[ind] += lnLContrib; 
                        partQty[ind] += lnLContrib;
                    }
                }
            }           
    }
    if(hasIonTerm){
        if(isPXiGrid){
            len_t np2_store = 1 + np2 - this->np2; // account for the +1 on p2 flux grid
            for(len_t i = 0; i<np1; i++){
                for(len_t ir = 0; ir<nr; ir++){
                    partContrib = preFactor[i]*lnLei[ir][pind];
                    for(len_t iz=0; iz<nZ; iz++)
                        for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                            indZ = ionIndex[iz][Z0]; 
                            real_t DpartContrib = ionDensities[ir][indZ] * preFactor[i] * lnLEI_partialNi[ir][indZ];
                            len_t Zfact;
                            if(isNonScreened)
                                Zfact = Zs[iz]*Zs[iz];
                            else 
                                Zfact = Z0*Z0;
                            real_t lnLContrib = Zfact*ionTerm[indZ*np1*np2_store+i]*DpartContrib;
                            real_t tmpQty = Zfact*ionTerm[indZ*np1*np2_store+i]*partContrib + lnLContrib;
                            for(len_t j = 0; j<np2; j++){
                                pind = np1*j+i;
                                len_t ind = (indZ*nr+ir)*np1*np2 + pind;
                                ionLnLContrib[ind] += lnLContrib;
                                partQty[ind] += tmpQty;
                            }
                    }
                }
            }
        } else {
            for(len_t i = 0; i<np1; i++)
                for(len_t j = 0; j<np2; j++){
                    pind = np1*j+i;
                    for(len_t ir = 0; ir<nr; ir++){
                        partContrib = preFactor[pind]*lnLei[ir][pind];
                        for(len_t iz=0; iz<nZ; iz++)
                            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                                indZ = ionIndex[iz][Z0]; 
                                len_t ind = (indZ*nr+ir)*np1*np2 + pind;
                                real_t DpartContrib = ionDensities[ir][indZ] * preFactor[pind] * lnLEI_partialNi[ir][indZ];
                                len_t Zfact;
                                if(isNonScreened)
                                    Zfact = Zs[iz]*Zs[iz];
                                else 
                                    Zfact = Z0*Z0;
                                real_t lnLContrib = Zfact*ionTerm[indZ*np1*np2+pind]*DpartContrib;
                                ionLnLContrib[ind] += lnLContrib;
                                partQty[ind] += Zfact*ionTerm[indZ*np1*np2+pind]*partContrib + lnLContrib;
                            }
                    }
                }
        }
    }
    if(isBrems){
        len_t np2_store;
        if(isPXiGrid)
            np2_store = 1 + np2 - this->np2; // account for the +1 on p2 flux grid
        else 
            np2_store = np2;
        for(len_t i = 0; i<np1; i++)
            for(len_t j = 0; j<np2; j++){
                pind = np1*j+i;
                if(isPXiGrid) 
                    pindStore = i;
                else 
                    pindStore = pind;
                
                for(len_t ir = 0; ir<nr; ir++)
                    for(len_t indZ=0; indZ<nzs; indZ++)
                        partQty[(indZ*nr+ir)*np1*np2 + pind] += bremsTerm[indZ*np1*np2_store+pindStore];
            }       
    } if(isNonScreened){
        for(len_t i = 0; i<np1; i++)
            for(len_t j = 0; j<np2; j++){
                pind = np1*j+i;
                if(isPXiGrid)
                    pindStore = i;
                else
                    pindStore = pind;
                
                for(len_t ir = 0; ir<nr; ir++){
                    electronTerm = nColdTerm[ir][pindStore]*preFactor[pindStore]*lnLee[ir][pind];
                    for(len_t iz=0; iz<nZ; iz++)
                        for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                            indZ = ionIndex[iz][Z0]; 
                            partQty[(indZ*nr+ir)*np1*np2 + pind] += (Zs[iz]-Z0)*electronTerm;
                        }
                }
            }
    } else if(isPartiallyScreened){
        if(isPXiGrid){
            len_t np2_store = 1 + np2 - this->np2; // account for the +1 on p2 flux grid
            for(len_t i = 0; i<np1; i++)
                for(len_t indZ=0; indZ<nzs; indZ++)
                    for(len_t ir = 0; ir<nr; ir++){
                        real_t tmpQty = preFactor[i]*screenedTerm[indZ*np1*np2_store + i];
                        for(len_t j = 0; j<np2; j++)
                            partQty[indZ*nr*np1*np2 + ir*np1*np2 + np1*j+i] += tmpQty;
                    }
        } else {
            for(len_t j = 0; j<np2; j++)
                for(len_t i = 0; i<np1; i++){
                    pind = np1*j+i;
                    for(len_t indZ=0; indZ<nzs; indZ++)
                        for(len_t ir = 0; ir<nr; ir++)
                            partQty[indZ*nr*np1*np2 + ir*np1*np2 + pind] += preFactor[pind]*screenedTerm[indZ*np1*np2 + pind];
                }
            }
        }
    
    for(len_t ir=0; ir<nr; ir++){
        delete [] lnLEE_partialNi[ir];
        delete [] lnLEI_partialNi[ir];
    }
    delete [] lnLEE_partialNi;
    delete [] lnLEI_partialNi;
}


void CollisionFrequency::SetNColdPartialContribution(real_t **nColdTerm,real_t *preFactor, real_t *const* lnLee, len_t nr, len_t np1, len_t np2, real_t *&partQty){
    if(partQty==nullptr)
        partQty = new real_t[np1*np2*nr];
    for(len_t it=0; it < np1*np2*nr; it++)
        partQty[it] = 0;

    len_t pind, pindStore;
    for(len_t i = 0; i<np1; i++)
        for(len_t j = 0; j<np2; j++){
            pind = np1*j+i;
            if(isPXiGrid)
                pindStore = i;
            else
                pindStore = pind;
            
            // TODO: Possible optimization: if(isPXiGrid), calculate RHS outside the j loop
            for(len_t ir = 0; ir<nr; ir++)
                partQty[np1*np2*ir + pind] = nColdTerm[ir][pindStore]*preFactor[pindStore]*lnLee[ir][pind];
        }
}


/**
 * Set partial derivative of quantity with respect to T_cold.
 * For now using a placeholder method where, if FULL operator,
 * assume a simple T^-1.5 dependence of the coefficient.
 */
void CollisionFrequency::SetTColdPartialContribution(real_t **nColdTerm, real_t *ionTerm, real_t *preFactor, real_t *const* lnLee,  const real_t *pIn, len_t nr, len_t np1, len_t np2, real_t *&TColdPartialContribution){
    if(TColdPartialContribution==nullptr)
        TColdPartialContribution = new real_t[np1*np2*nr];    
    for(len_t it=0; it < np1*np2*nr; it++)
        TColdPartialContribution[it] = 0;

    len_t pind;
    
    const real_t *ncold = unknowns->GetUnknownData(id_ncold);
    if(isPXiGrid)
        for(len_t i=0;i<np1;i++)
            for(len_t ir=0; ir<nr; ir++){
                real_t DDTElectronTerm = evaluateDDTElectronTermAtP(ir,pIn[i],collQtySettings->collfreq_mode);
                real_t dLnL = lnLambdaEE->evaluatePartialAtP(ir,pIn[i],id_Tcold,0);
                for(len_t j=0;j<np2;j++){
                    pind = np1*j+i;
                    TColdPartialContribution[np1*np2*ir + pind] = ncold[ir] * preFactor[i] *
                        (lnLee[ir][pind]*DDTElectronTerm + dLnL * nColdTerm[ir][i]);
                }
            }
    else
        for(len_t pind=0;pind<np1*np2;pind++)
            for(len_t ir=0; ir<nr; ir++){
                real_t dLnL = lnLambdaEE->evaluatePartialAtP(ir,pIn[pind],id_Tcold,0);
                TColdPartialContribution[np1*np2*ir + pind] = ncold[ir]*preFactor[pind] * 
                    (lnLee[ir][pind]*evaluateDDTElectronTermAtP(ir,pIn[pind],collQtySettings->collfreq_mode) + dLnL * nColdTerm[ir][pind]);
            }
    len_t indZ, Zfact;
    if(hasIonTerm){
        if(isPXiGrid){
            len_t np2_store = 1 + np2 - this->np2; // account for the +1 on p2 flux grid
            for(len_t i=0;i<np1;i++)
                for(len_t ir=0; ir<nr; ir++){
                    real_t lnLEI_partialT = lnLambdaEI->evaluatePartialAtP(ir, pIn[i], id_Tcold, 0);
                    for(len_t iz = 0; iz<nZ; iz++)
                        for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                            indZ = ionIndex[iz][Z0];
                            if(isNonScreened)
                                Zfact = Zs[iz]*Zs[iz];
                            else
                                Zfact = Z0*Z0;
                            real_t PZFactor = Zfact * preFactor[i] * ionTerm[indZ*np1*np2_store+i];
                            real_t TCold_tmp = PZFactor * lnLEI_partialT * ionDensities[ir][indZ];
                            for(len_t j=0;j<np2;j++)
                                TColdPartialContribution[np1*np2*ir + np1*j+i] += TCold_tmp;
                        }
                }
        } else
            for(len_t i=0;i<np1;i++)
                for(len_t j=0;j<np2;j++){
                    pind = np1*j+i;
                    for(len_t iz = 0; iz<nZ; iz++)
                        for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                            indZ = ionIndex[iz][Z0];
                            if(isNonScreened)
                                Zfact = Zs[iz]*Zs[iz];
                            else
                                Zfact = Z0*Z0;
                            real_t PZFactor = Zfact * preFactor[pind] * ionTerm[indZ*np1*np2+pind];
                            for(len_t ir=0; ir<nr; ir++)
                                TColdPartialContribution[np1*np2*ir + pind] += PZFactor * lnLambdaEI->evaluatePartialAtP(ir, pIn[pind], id_Tcold, 0) * ionDensities[ir][indZ];
                        }
                }
            }
}


/**
 * Sets the partial derivative of the frequency with respect to f_hot. 
 * (i.e. the distribution function on a hot-tail grid)
 */
void CollisionFrequency::SetNonlinearPartialContribution(CoulombLogarithm *lnLambda, real_t *&partQty){
    if(partQty==nullptr)
        partQty = new real_t[np1*(np1+1)*nr];
    for(len_t it=0; it < np1*(np1+1)*nr; it++)
        partQty[it] = 0;
    for(len_t i=0; i<np1+1; i++)
        for(len_t ir=0;ir<nr;ir++)
            for(len_t ip=0; ip<np1; ip++)
                partQty[ip*(np1+1)*nr + (np1+1)*ir + i] = lnLambda->GetLnLambdaT(ir)*nonlinearMat[i][ip];
}


/**
 * Allocates quantities which will be used in the calculation of the collision frequencies.
 */
void CollisionFrequency::AllocatePartialQuantities(){
    DeallocatePartialQuantities();
    InitializeGSLWorkspace();
    Zs = new real_t[nZ];
    ionIndex = new real_t*[nZ];
    ionDensities = new real_t*[nr];
    atomicParameter = new real_t[nzs];

    for(len_t iz=0;iz<nZ;iz++)
        ionIndex[iz] = new real_t[ionHandler->GetZ(iz)+1];
    for(len_t ir=0; ir<nr;ir++)
        ionDensities[ir] = new real_t[nzs];
    if(!buildOnlyF1F2){
        preFactor    = new real_t[np1*np2_store];
        preFactor_fr = new real_t[np1*np2_store];
        if(hasIonTerm){
            ionTerm = new real_t[nzs*np1*np2_store];
            ionTerm_fr = new real_t[nzs*np1*np2_store];
        }
        ionPartialContribution    = new real_t[nzs*nr*np1*np2];
        ionPartialContribution_fr = new real_t[nzs*(nr+1)*np1*np2];
        ionLnLambdaPartialContribution    = new real_t[nzs*nr*np1*np2];
        ionLnLambdaPartialContribution_fr = new real_t[nzs*(nr+1)*np1*np2];
        
        if(isPartiallyScreened){
            screenedTerm    = new real_t[nzs*np1*np2_store];
            screenedTerm_fr = new real_t[nzs*np1*np2_store];
        }
        if(isBrems){
            bremsTerm    = new real_t[nzs*np1*np2_store];
            bremsTerm_fr = new real_t[nzs*np1*np2_store];
        }
        nColdTerm = new real_t*[nr];
        nColdTerm_fr = new real_t*[nr+1];
        for(len_t ir=0;ir<nr;ir++)
            nColdTerm[ir] = new real_t[np1*np2_store];
        for(len_t ir=0;ir<nr+1;ir++)
            nColdTerm_fr[ir] = new real_t[np1*np2_store];
        nColdPartialContribution    = new real_t[nr*np1*np2];
        nColdPartialContribution_fr = new real_t[(nr+1)*np1*np2];
        TColdPartialContribution    = new real_t[nr*np1*np2];
        TColdPartialContribution_fr = new real_t[(nr+1)*np1*np2];
    }
    preFactor_f1 = new real_t[(np1+1)*np2_store];
    preFactor_f2 = new real_t[np1*(np2_store+1)];
    if(hasIonTerm){
        ionTerm_f1 = new real_t[nzs*(np1+1)*np2_store];
        ionTerm_f2 = new real_t[nzs*np1*(np2_store+1)];
    }
    ionPartialContribution_f1 = new real_t[nzs*nr*(np1+1)*np2];
    ionPartialContribution_f2 = new real_t[nzs*nr*np1*(np2+1)];
    ionLnLambdaPartialContribution_f1 = new real_t[nzs*nr*(np1+1)*np2];
    ionLnLambdaPartialContribution_f2 = new real_t[nzs*nr*np1*(np2+1)];
    if(isPartiallyScreened){
        screenedTerm_f1 = new real_t[nzs*(np1+1)*np2_store];
        screenedTerm_f2 = new real_t[nzs*np1*(np2_store+1)];
    }
    if(isBrems){
        bremsTerm_f1 = new real_t[nzs*(np1+1)*np2_store];
        bremsTerm_f2 = new real_t[nzs*np1*(np2_store+1)];
    }
    nColdTerm_f1 = new real_t*[nr];
    nColdTerm_f2 = new real_t*[nr];
    for(len_t ir=0;ir<nr;ir++){
        nColdTerm_f1[ir] = new real_t[(np1+1)*np2_store];
        nColdTerm_f2[ir] = new real_t[np1*(np2_store+1)];
    }
    nColdPartialContribution_f1 = new real_t[nr*(np1+1)*np2];
    nColdPartialContribution_f2 = new real_t[nr*np1*(np2+1)];
    TColdPartialContribution_f1 = new real_t[nr*(np1+1)*np2];
    TColdPartialContribution_f2 = new real_t[nr*np1*(np2+1)];

    if (isNonlinear){
        nonlinearMat = new real_t*[np1+1]; // multiply matrix by f lnLc to get p*nu_s on p flux grid
        for (len_t i = 0; i<np1+1; i++)
            nonlinearMat[i] = new real_t[np1];

        const real_t *p = mg->GetP1();
        trapzWeights = new real_t[np1];
        for (len_t i = 1; i<np1-1; i++)
            trapzWeights[i] = (p[i+1]-p[i-1])/2;

        fHotPartialContribution_f1 = new real_t[nr*np1*(np1+1)];    
    }
}


/**
 * Deallocator
 */
void CollisionFrequency::DeallocatePartialQuantities(){
    DeallocateGSL();
    if (Zs != nullptr){
        delete [] Zs;
        for(len_t iz=0; iz<nZ; iz++)
            delete [] ionIndex[iz];
        for(len_t ir=0; ir<nr;ir++)
            delete [] ionDensities[ir];
        
        delete [] ionIndex;
        delete [] ionDensities; 
    }
    if(preFactor!=nullptr){
        delete [] preFactor;
        delete [] preFactor_fr;
    }
    if(ionTerm!=nullptr){
        delete [] ionTerm;
        delete [] ionTerm_fr;
    }   
    if(preFactor_f1 != nullptr){
        delete [] preFactor_f1;
        delete [] preFactor_f2;
    }
    if(ionTerm!=nullptr){
        delete [] ionTerm_f1;
        delete [] ionTerm_f2;
    }
    if(nColdTerm != nullptr){
        for(len_t ir=0;ir<nr;ir++)
            delete [] nColdTerm[ir];
        for(len_t ir=0;ir<nr+1;ir++)
            delete [] nColdTerm_fr[ir];
        delete [] nColdTerm;
        delete [] nColdTerm_fr;
    }
    if (screenedTerm != nullptr){
        delete [] screenedTerm;
        delete [] screenedTerm_fr;            
    }
    if (screenedTerm_f1 != nullptr){
        delete [] screenedTerm_f1;
        delete [] screenedTerm_f2;
    }
    if (nColdPartialContribution != nullptr){
        delete [] nColdPartialContribution;
        delete [] nColdPartialContribution_fr;
    }
    if (nColdPartialContribution_f1 != nullptr){
        delete [] nColdPartialContribution_f1;
        delete [] nColdPartialContribution_f2;
    }
    if (TColdPartialContribution != nullptr){
        delete [] TColdPartialContribution;
        delete [] TColdPartialContribution_fr;
    }
    if (TColdPartialContribution_f1 != nullptr){
        delete [] TColdPartialContribution_f1;
        delete [] TColdPartialContribution_f2;
    }
    if (ionPartialContribution != nullptr){
        delete [] ionPartialContribution;
        delete [] ionPartialContribution_fr;
        delete [] ionLnLambdaPartialContribution;
        delete [] ionLnLambdaPartialContribution_fr;
    }
    if (ionPartialContribution_f1 != nullptr){
        delete [] ionPartialContribution_f1;
        delete [] ionPartialContribution_f2;
        delete [] ionLnLambdaPartialContribution_f1;
        delete [] ionLnLambdaPartialContribution_f2;
    }
    if(nColdTerm_f1 != nullptr){
        for(len_t ir=0;ir<nr;ir++){
            delete [] nColdTerm_f1[ir];
            delete [] nColdTerm_f2[ir];
        }
        delete [] nColdTerm_f1;
        delete [] nColdTerm_f2;
    }
    if(atomicParameter != nullptr)
        delete [] atomicParameter;
    if(nonlinearMat != nullptr){
        for(len_t i = 0; i<np1+1;i++)
            delete [] nonlinearMat[i];
        delete [] nonlinearMat;
        delete [] trapzWeights;
        delete [] fHotPartialContribution_f1;
    }
}


/**
 * Initializes a GSL workspace for evaluation of full e-e test particle operator
 */
void CollisionFrequency::InitializeGSLWorkspace(){
    DeallocateGSL();
    gsl_ad_w = gsl_integration_workspace_alloc(1000); 
}


/**
 * Deallocator
 */
void CollisionFrequency::DeallocateGSL(){
    if (this->gsl_ad_w != nullptr)
        gsl_integration_workspace_free(gsl_ad_w);
}


/**
 * Evaluates the Jacobian with respect to unknown derivId of the  
 * collision frequency at radial grid point ir and momentum p.
 */
real_t CollisionFrequency::evaluatePartialAtP(len_t ir, real_t p, len_t derivId, len_t n,struct collqty_settings *inSettings){     
    // Return 0 for all other derivId but ncold, ni or Tcold (when collfreq mode is FULL)
    if( ! ( derivId == id_ncold || derivId == id_ni || derivId == id_Tcold  ) )
        return 0;

    bool isPartiallyScreened = (inSettings->collfreq_type == OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED);
    bool isNonScreened = (inSettings->collfreq_type == OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED);
    bool isBrems = (inSettings->bremsstrahlung_mode != OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT);


    real_t preFact = evaluatePreFactorAtP(p,inSettings->collfreq_mode); 
    real_t lnLee = lnLambdaEE->evaluateAtP(ir,p,inSettings);
    real_t lnLei = lnLambdaEI->evaluateAtP(ir,p,inSettings);
    
    real_t electronTerm = evaluateElectronTermAtP(ir,p,inSettings->collfreq_mode);

    // if ncold, this is the jacobian
    if(derivId == id_ncold)
        return preFact*lnLee *electronTerm;
    // else, for ions or Tcold, we move on...

    real_t dLnLee = lnLambdaEE->evaluatePartialAtP(ir,p,derivId,n,inSettings);
    real_t dLnLei = lnLambdaEI->evaluatePartialAtP(ir,p,derivId,n,inSettings);
    
    real_t ntarget = unknowns->GetUnknownData(id_ncold)[ir];
    if (isNonScreened)
        ntarget += ionHandler->GetBoundElectronDensity(ir);
    // evaluate and return T_cold expression
    if(derivId == id_Tcold){
        real_t DDTelectronTerm = lnLee*evaluateDDTElectronTermAtP(ir,p,inSettings->collfreq_mode) 
                                + dLnLee*electronTerm;
        real_t electronContribution = preFact * ntarget * DDTelectronTerm; 
        real_t ionContribution = 0;
        if(hasIonTerm)
            for(len_t iz = 0; iz<nZ; iz++)
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                    len_t Zfact = Z0*Z0;
                    len_t indZ = ionIndex[iz][Z0];
                    if(isNonScreened)
                        Zfact = Zs[iz]*Zs[iz];
                    ionContribution += preFact * dLnLei * Zfact * evaluateIonTermAtP(iz,Z0,p) * ionDensities[ir][indZ];
                }
        return electronContribution + ionContribution;
    }

    // else treat n_i case
    // set iz and Z0 corresponding to input nMultiple "n"
    len_t iz_in, Z0_in;
    ionHandler->GetIonIndices(n, iz_in, Z0_in);
    len_t Zs_in = Zs[iz_in];

    real_t collFreq = dLnLee*electronTerm*ntarget;
    if(isNonScreened)
        collFreq += (Zs_in-Z0_in) * lnLee *electronTerm;

    // Add ion contribution; SlowingDownFrequency doesn't have one and will skip this step
    if(hasIonTerm){
        len_t Zfact;
        if(isNonScreened)
            Zfact = Zs_in*Zs_in;
        else 
            Zfact = Z0_in*Z0_in;
        collFreq += (lnLei + dLnLei*ionDensities[ir][n]) * Zfact * evaluateIonTermAtP(iz_in,Z0_in,p);
    }
    // Add screening contribution
    if(isPartiallyScreened)
        collFreq += evaluateScreenedTermAtP(iz_in,Z0_in,p,inSettings->collfreq_mode);

    collFreq *= preFact;

    // Add Bremsstrahlung contribution
    if(isBrems)
        collFreq += evaluateBremsstrahlungTermAtP(iz_in,Z0_in,p,inSettings->bremsstrahlung_mode,inSettings->collfreq_type);

    return collFreq;
}
