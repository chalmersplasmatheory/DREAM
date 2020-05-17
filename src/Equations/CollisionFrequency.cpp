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

using namespace DREAM;

CollisionFrequency::CollisionFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                CoulombLogarithm *lnLee, CoulombLogarithm *lnLei,
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantityHandler::collqtyhand_settings *cqset)
                : CollisionQuantity(g,u,ih,mgtype,cqset) {
    lnLambdaEE = lnLee;
    lnLambdaEI = lnLei;
    
}

/**
 * Evaluates the collision frequency at radial grid point ir and momentum p,
 * neglecting any contribution from the nonlinear collision operator
 */
real_t CollisionFrequency::evaluateAtP(len_t ir, real_t p){
    real_t *ncold = unknowns->GetUnknownData(id_ncold);
    real_t ntarget = ncold[ir];
    if (isNonScreened)
        ntarget += nbound[ir];

    real_t preFact = evaluatePreFactorAtP(p); 
    real_t lnLee = lnLambdaEE->evaluateAtP(ir,p);
    real_t lnLei = lnLambdaEI->evaluateAtP(ir,p);
    
    // Add electron contribution to collision frequency
    real_t collQty = lnLee * evaluateElectronTermAtP(ir,p) * ntarget;

    len_t ind;
    // Add ion contribution 
    if(hasIonTerm){
        if(isNonScreened){
            for(len_t iz = 0; iz<nZ; iz++){
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                    ind = ionIndex[iz][Z0];
                    collQty += lnLei * Zs[iz]*Zs[iz] * evaluateIonTermAtP(iz,Z0,p) * ionDensities[ir][ind];
                }
            }

        } else {
            for(len_t iz = 0; iz<nZ; iz++){
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                    ind = ionIndex[iz][Z0];
                    collQty += lnLei * Z0*Z0 * evaluateIonTermAtP(iz,Z0,p) * ionDensities[ir][ind];
                }
            }
        }
    }
    // Add screening contribution
    if(isPartiallyScreened){
        for(len_t iz = 0; iz<nZ; iz++){
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                ind = ionIndex[iz][Z0];
                collQty +=  evaluateScreenedTermAtP(iz,Z0,p) * ionDensities[ir][ind];
            }
        }
    }
    collQty *= preFact;
    return collQty;
}

// Rebuilds partial terms which depend on unknown quantities (density and temperature).
void CollisionFrequency::RebuildPlasmaDependentTerms(){
    nbound = ionHandler->evaluateBoundElectronDensityFromQuasiNeutrality(nbound);
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
        setNColdTerm(nColdTerm,mg->GetP(),nr,np1,np2_store);
        setNColdTerm(nColdTerm_fr,mg->GetP(),nr+1,np1,np2_store);
    }
    setNColdTerm(nColdTerm_f1,mg->GetP_f1(),nr,np1+1,np2_store);
    setNColdTerm(nColdTerm_f2,mg->GetP_f2(),nr,np1,np2_store+1);
}


/**
 * Rebuilds partial contributions that only depend on the grid. If using P-Xi grid, 
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
        }
    }
    len_t ind;
    for(len_t iZ = 0; iZ<nZ; iZ++){
        for(len_t Z0=0; Z0<=Zs[iZ]; Z0++){
            ind = ionIndex[iZ][Z0];
            atomicParameter[ind] = GetAtomicParameter(iZ,Z0);
        }
    }

    if (!buildOnlyF1F2){
        setPreFactor(preFactor,mg->GetP(),np1,np2_store);
        setPreFactor(preFactor_fr,mg->GetP(),np1,np2_store);
        setIonTerm(ionTerm,mg->GetP(),np1,np2_store);
        setIonTerm(ionTerm_fr,mg->GetP(),np1,np2_store);
        if(isPartiallyScreened){
            setScreenedTerm(screenedTerm,mg->GetP(),np1,np2_store);
            setScreenedTerm(screenedTerm_fr,mg->GetP(),np1,np2_store);
        }
    }
    setPreFactor(preFactor_f1,mg->GetP_f1(),np1+1,np2_store);
    setPreFactor(preFactor_f2,mg->GetP_f2(),np1,np2_store+1);
    setIonTerm(ionTerm_f1,mg->GetP_f1(),np1+1,np2_store);
    setIonTerm(ionTerm_f2,mg->GetP_f2(),np1,np2_store+1);
    if(isPartiallyScreened){
        setScreenedTerm(screenedTerm_f1,mg->GetP_f1(),np1+1,np2_store);
        setScreenedTerm(screenedTerm_f2,mg->GetP_f2(),np1,np2_store+1);
    }
    if(isNonlinear)
        calculateIsotropicNonlinearOperatorMatrix();
}

/**
 * Calculates and stores the momentum-dependent prefactor to the collision frequencies.
 */
void CollisionFrequency::setPreFactor(real_t *&preFactor, const real_t *pIn, len_t np1, len_t np2){
    real_t p, PF;
    len_t ind;
    for(len_t i = 0; i<np1; i++){
        for (len_t j = 0; j<np2; j++){
            ind = np1*j+i;
            p = pIn[ind];
            if(p==0)
                PF = ReallyLargeNumber; 
            else
                PF = evaluatePreFactorAtP(p);
            preFactor[ind] = PF;
        }
    }
}

/**
 * Puts together all the partial contributions to the collision frequency to get the full thing.
 */
void CollisionFrequency::AssembleQuantity(real_t **&collisionQuantity, len_t nr, len_t np1, len_t np2, len_t fluxGridType){
    real_t *nColdContribution = new real_t[nr*np1*np2];
    real_t *niContribution = new real_t[nzs*nr*np1*np2];
    real_t collQty;
    real_t *ncold = unknowns->GetUnknownData(id_ncold);
    const len_t *Zs = ionHandler->GetZs();
    nColdContribution = GetUnknownPartialContribution(id_ncold,fluxGridType,nColdContribution);
    niContribution    = GetUnknownPartialContribution(id_ni,   fluxGridType,niContribution);

    len_t indZ;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t j=0; j<np2; j++){
            for(len_t i=0; i<np1; i++){
                collQty = ncold[ir]*nColdContribution[np1*np2*ir + np1*j + i];
                for(len_t iz = 0; iz<nZ; iz++){
                    for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                        indZ = ionIndex[iz][Z0];
                        collQty += ionDensities[ir][indZ]*niContribution[indZ*nr*np1*np2 + np1*np2*ir + np1*j + i];
                    }
                }
                collisionQuantity[ir][j*np1+i] = collQty; 
            }
        }
    }
    delete [] nColdContribution;
    delete [] niContribution;
}

/**
 *  Calculation of the partial contribution to the collision frequency from the unknown quantity
 * with ID id_unknown. Returns the partial derivative of the term with respect to that quantity 
 * (ignoring variations with lnLambda). See AssembleQuantity(...) or addNonlinearContribution() 
 * for how it is used.
 */

real_t *CollisionFrequency::GetUnknownPartialContribution(len_t id_unknown, len_t fluxGridMode, real_t *&partQty){
    if(id_unknown == id_ncold)
        GetNColdPartialContribution(fluxGridMode,partQty);
    else if(id_unknown == id_ni)
        GetNiPartialContribution(fluxGridMode,partQty);
    else if(id_unknown == id_fhot){
        if(!( (fluxGridMode==2)&&(np2=1)&&(isPXiGrid) ) )
            throw FVM::FVMException("Nonlinear contribution to collision frequencies is only implemented for hot-tails, with p-xi grid and np2=1 and evaluated on the p flux grid.");
        GetNonlinearPartialContribution(lnLambdaEE->GetLnLambdaC(),partQty);

    } else
        throw FVM::FVMException("Invalid id_unknown: %s does not contribute to the collision frequencies",unknowns->GetUnknown(id_unknown)->GetName());
    return partQty;
}


/** Adds the non-linear contribution to the collision frequency. For now, only supports 
 * hot-tail grids where np2=1 and using a pxi-grid, and only updates the p flux grid 
 * component.
 */
void CollisionFrequency::AddNonlinearContribution(){
    real_t *fHot = unknowns->GetUnknownData(id_fhot);
    real_t *fHotContribution = new real_t[nr*np1*(np1+1)];
    GetNonlinearPartialContribution(lnLambdaEE->GetLnLambdaC(),fHotContribution);

    for (len_t ir=0;ir<nr;ir++)
        for(len_t i=0; i<np1+1; i++)
            for(len_t ip=0; ip<np1; ip++)
                collisionQuantity_f1[ir][i] += fHotContribution[ip*(np1+1)*nr + ir*(np1+1) + i] * fHot[np1*ir+ip];
}


/**
 * Calculates and stores the ion contribution to the collision frequency.
 */
void CollisionFrequency::setIonTerm(real_t *&ionTerm, const real_t *pIn, len_t np1, len_t np2){
    if(!hasIonTerm)
        return;

    real_t p;
    len_t ind, pind;
    for(len_t i = 0; i<np1; i++){
        for (len_t j = 0; j<np2; j++){
            pind = np1*j+i;
            p = pIn[pind];
            for(len_t iz = 0; iz<nZ; iz++){
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                    ind = ionIndex[iz][Z0];
                    ionTerm[ind*np1*np2 + pind] = evaluateIonTermAtP(iz,Z0,p);
                }
            }
        }
    }
}

/**
 * Calculates and stores the partially-screened contribution to the collision frequency.
 */
void CollisionFrequency::setScreenedTerm(real_t *&screenedTerm, const real_t *pIn, len_t np1, len_t np2){

    real_t p;
    len_t ind, pind;
    for(len_t i = 0; i<np1; i++){
        for (len_t j = 0; j<np2; j++){
            pind = np1*j+i;
            p = pIn[pind];
            for(len_t iz = 0; iz<nZ; iz++){
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                    ind = ionIndex[iz][Z0];
                    screenedTerm[ind*np1*np2 + pind] = evaluateScreenedTermAtP(iz,Z0,p);
                }
            }
        }
    }
}

/**
 * Calculates and stores the free-electron contribution to the collision frequency.
 */
void CollisionFrequency::setNColdTerm(real_t **&nColdTerm, const real_t *pIn, len_t nr, len_t np1, len_t np2){
    real_t p;
    len_t pind;
    // Depending on setting, set nu_s to superthermal or full formula (with maxwellian)
    if (collQtySettings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL) {
        for(len_t i=0;i<np1;i++){
            for(len_t j=0;j<np2;j++){
                pind = np1*j+i;
                p = pIn[pind];
                for(len_t ir=0; ir<nr; ir++)
                    nColdTerm[ir][pind] = evaluateElectronTermAtP(ir,p);
            }
        }
    } 
}




/**
 * The following methods are helper functions for the evaluation of
 * the frequencies appearing in the relativistic test-particle operator. 
 */
real_t CollisionFrequency::psi0Integrand(real_t x, void *params){
    real_t gamma = *(real_t *) params;
    return 1/sqrt( (x+gamma)*(x+gamma)-1 );
} 
real_t CollisionFrequency::psi1Integrand(real_t x, void *params){
    real_t gamma = *(real_t *) params;
    return (x+gamma)/sqrt((x+gamma)*(x+gamma)-1); // integrated with weight w(x) = exp(-(x-gamma)/Theta) 
} 
/** 
 * Evaluates integral appearing in relativistic test-particle operator
 * Psi0 = int_0^p exp( -(sqrt(1+s^2)-1)/Theta) / sqrt(1+s^2) ds;
 */
real_t CollisionFrequency::evaluatePsi0(len_t ir, real_t p) {
    real_t gamma = sqrt(1+p*p);
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    
    gsl_function F;
    F.function = &(CollisionFrequency::psi0Integrand); 
    F.params = &gamma;
    real_t psi0int; 
    gsl_integration_fixed(&F, &psi0int, gsl_w[ir]);

    real_t Theta = T_cold[ir] / Constants::mc2inEV;
    return evaluateExp1OverThetaK(Theta,0) - exp( -(gamma-1)/Theta ) * psi0int;
}
real_t CollisionFrequency::evaluatePsi1(len_t ir, real_t p) {
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    real_t gamma = sqrt(1+p*p);
    gsl_function F;
    F.function = &(CollisionFrequency::psi1Integrand); 
    F.params = &gamma;
    real_t psi1int; 
    gsl_integration_fixed(&F, &psi1int, gsl_w[ir]);

    real_t Theta = T_cold[ir] / Constants::mc2inEV;
    return evaluateExp1OverThetaK(Theta,1) - exp( -(gamma-1)/Theta ) * psi1int;
}


// evaluates e^x K_n(x), with K_n the (exponentially decreasing) modified bessel function.
real_t CollisionFrequency::evaluateExp1OverThetaK(real_t Theta, real_t n) {
    real_t ThetaThreshold = 0.002;
    /**
     * Since cyl_bessel_k ~ exp(-1/Theta), for small Theta you get precision issues.
     * Instead using the asymptotic expansion for bessel_k for small Theta.
     */
    if (Theta > ThetaThreshold)
        return exp(1/Theta)*std::cyl_bessel_k(n,1/Theta);
    else {
        real_t n2 = n*n;
        return sqrt(M_PI*Theta/2)*(1 + (4*n2-1)/8 * Theta + (4*n2-1)*(4*n2-9)*Theta*Theta/128 + (4*n2-1)*(4*n2-9)*(4*n2-25)*Theta*Theta*Theta/3072);
    }
}


/**
 * Returns the partial derivative of the frequency with respect to n_cold.
 */
void CollisionFrequency::GetNColdPartialContribution(len_t fluxGridMode, real_t *&partQty){
    if(fluxGridMode==0){
        GetNColdPartialContribution(nColdTerm,preFactor,lnLambdaEE->GetValue(),nr,np1,np2,partQty);
    } else if(fluxGridMode==1){
        GetNColdPartialContribution(nColdTerm_fr,preFactor_fr,lnLambdaEE->GetValue_fr(),nr+1,np1,np2,partQty);
    } else if(fluxGridMode==2){
        GetNColdPartialContribution(nColdTerm_f1,preFactor_f1,lnLambdaEE->GetValue_f1(),nr,np1+1,np2,partQty);
    } else if(fluxGridMode==3){
        GetNColdPartialContribution(nColdTerm_f2,preFactor_f2,lnLambdaEE->GetValue_f2(),nr,np1,np2+1,partQty);
    } else
        throw FVM::FVMException("Invalid fluxGridMode.");
}
/**
 * Returns the partial derivative of the frequency with respect to n_i.
 */
void CollisionFrequency::GetNiPartialContribution(len_t fluxGridMode, real_t *&partQty){
    if(fluxGridMode==0){
        GetNiPartialContribution(nColdTerm,ionTerm, screenedTerm,preFactor,lnLambdaEE->GetValue(),lnLambdaEI->GetValue(),nr,np1,np2,partQty);
    } else if(fluxGridMode==1){
        GetNiPartialContribution(nColdTerm_fr,ionTerm_fr,screenedTerm_fr,preFactor_fr,lnLambdaEE->GetValue_fr(),lnLambdaEI->GetValue_fr(),nr+1,np1,np2,partQty);
    } else if(fluxGridMode==2){
        GetNiPartialContribution(nColdTerm_f1,ionTerm_f1,screenedTerm_f1,preFactor_f1,lnLambdaEE->GetValue_f1(),lnLambdaEI->GetValue_f1(),nr,np1+1,np2,partQty);
    } else if(fluxGridMode==3){
        GetNiPartialContribution(nColdTerm_f2,ionTerm_f2,screenedTerm_f2,preFactor_f2,lnLambdaEE->GetValue_f2(),lnLambdaEI->GetValue_f2(),nr,np1,np2+1,partQty);
    } else
        throw FVM::FVMException("Invalid fluxGridMode.");
}


void CollisionFrequency::GetNiPartialContribution(real_t **nColdTerm, real_t *ionTerm, real_t *screenedTerm, real_t *preFactor, real_t *const* lnLee,  real_t *const* lnLei, len_t nr, len_t np1, len_t np2, real_t *&partQty){
    if(partQty==nullptr){
        partQty = new real_t[nzs*np1*np2*nr];
    }
    // TODO: look over code to ensure that you have something freshly allocated so that this isn't needed    
    for(len_t it = 0; it<nzs*np1*np2*nr; it++){
        partQty[it] = 0;
    }
    len_t pind, pindStore, indZ;

    if(hasIonTerm){
        for(len_t i = 0; i<np1; i++){
            for(len_t j = 0; j<np2; j++){
                pind = np1*j+i;
                if(isPXiGrid) {
                    pindStore = i;
                } else {
                    pindStore = pind;
                }
                
                for(len_t ir = 0; ir<nr; ir++){
                    for(len_t iz=0; iz<nZ; iz++){
                        for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                            indZ = ionIndex[iz][Z0]; 
                            partQty[indZ*nr*np1*np2 + ir*np1*np2 + pind] = Z0*Z0*ionTerm[indZ*np1*np2_store+pindStore]*preFactor[pindStore]*lnLei[ir][pind];
                        }
                    }
                }
            }
        }
    }

    if(isNonScreened){
        real_t electronTerm;
        for(len_t i = 0; i<np1; i++){
            for(len_t j = 0; j<np2; j++){
                pind = np1*j+i;
                if(isPXiGrid)
                    pindStore = i;
                else
                    pindStore = pind;
                
                for(len_t ir = 0; ir<nr; ir++){
                    electronTerm = nColdTerm[ir][pindStore]*preFactor[pindStore]*lnLee[ir][pind];
                    for(len_t iz=0; iz<nZ; iz++){
                        for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                            indZ = ionIndex[iz][Z0]; 
                            partQty[indZ*nr*np1*np2 + ir*np1*np2 + pind] += (Zs[iz]-Z0)*electronTerm;
                        }
                    }
                    if(hasIonTerm){
                        for(len_t iz=0; iz<nZ; iz++){
                            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                                indZ = ionIndex[iz][Z0]; 
                                partQty[indZ*nr*np1*np2 + ir*np1*np2 + pind] += (Zs[iz]*Zs[iz]-Z0*Z0)*ionTerm[indZ*np1*np2_store+pindStore]*preFactor[pindStore]*lnLei[ir][pind];
                            }
                        }
                    }
                }
            }
        }
    } else if(isPartiallyScreened){
        for(len_t j = 0; j<np2; j++){
            for(len_t i = 0; i<np1; i++){
                pind = np1*j+i;
                if(isPXiGrid)
                    pindStore = i;
                 else 
                    pindStore = pind;
                
                
                for(len_t iz=0; iz<nZ; iz++){
                    for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                        indZ = ionIndex[iz][Z0]; 
                        for(len_t ir = 0; ir<nr; ir++){
                            partQty[indZ*nr*np1*np2 + ir*np1*np2 + pind] += preFactor[pindStore]*screenedTerm[indZ*np1*np2_store + pindStore];
                        }
                    }
                }
            }
        }

    }
}

void CollisionFrequency::GetNColdPartialContribution(real_t **nColdTerm,real_t *preFactor, real_t *const* lnLee, len_t nr, len_t np1, len_t np2, real_t *&partQty){
    if(partQty==nullptr){
        partQty = new real_t[np1*np2*nr];
    }
    for(len_t it=0; it < np1*np2*nr; it++){
        partQty[it] = 0;
    }
    len_t pind, pindStore;
    for(len_t i = 0; i<np1; i++){
        for(len_t j = 0; j<np2; j++){
            pind = np1*j+i;
            if(isPXiGrid)
                pindStore = i;
            else
                pindStore = pind;
            
            for(len_t ir = 0; ir<nr; ir++){
                partQty[np1*np2*ir + pind] = nColdTerm[ir][pindStore]*preFactor[pindStore]*lnLee[ir][pind];
            }
        }
    }
}



/**
 * Returns the partial derivative of the frequency with respect to f_hot. 
 * (i.e. the distribution function on a hot-tail grid)
 */
void CollisionFrequency::GetNonlinearPartialContribution(const real_t* lnLc, real_t *&partQty){
    if(partQty==nullptr){
        partQty = new real_t[np1*(np1+1)*nr];
    }

    for(len_t it=0; it < np1*(np1+1)*nr; it++){
        partQty[it] = 0;
    }

    for(len_t i=0; i<np1+1; i++)
        for(len_t ir=0;ir<nr;ir++)
            for(len_t ip=0; ip<np1; ip++)
                partQty[ip*(np1+1)*nr + (np1+1)*ir + i] = lnLc[ir]*nonlinearMat[i][ip];
}






/**
 * Allocates quantities which will be used in the calculation of the collision frequencies.
 */
void CollisionFrequency::AllocatePartialQuantities(){
    
    DeallocatePartialQuantities();
    InitializeGSLWorkspace();
    nbound = new real_t[nr];
    Zs = new real_t[nZ];
    ionIndex = new real_t*[nZ];
    ionDensities = new real_t*[nr];
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
    }
    preFactor_f1 = new real_t[(np1+1)*np2_store];
    preFactor_f2 = new real_t[np1*(np2_store+1)];
    if(hasIonTerm){
        ionTerm_f1 = new real_t[nzs*(np1+1)*np2_store];
        ionTerm_f2 = new real_t[nzs*np1*(np2_store+1)];
    }
    if(!buildOnlyF1F2){
        if(isPartiallyScreened){
            screenedTerm    = new real_t[nzs*np1*np2_store];
            screenedTerm_fr = new real_t[nzs*np1*np2_store];
        }
        nColdTerm = new real_t*[nr];
        nColdTerm_fr = new real_t*[nr+1];
        for(len_t ir=0;ir<nr;ir++)
            nColdTerm[ir] = new real_t[np1*np2_store];
        for(len_t ir=0;ir<nr+1;ir++)
            nColdTerm_fr[ir] = new real_t[np1*np2_store];
    }
        
        
    if(isPartiallyScreened){
        screenedTerm_f1 = new real_t[nzs*(np1+1)*np2_store];
        screenedTerm_f2 = new real_t[nzs*np1*(np2_store+1)];
    }

    nColdTerm_f1 = new real_t*[nr];
    nColdTerm_f2 = new real_t*[nr];
    for(len_t ir=0;ir<nr;ir++){
        nColdTerm_f1[ir] = new real_t[(np1+1)*np2_store];
        nColdTerm_f2[ir] = new real_t[np1*(np2_store+1)];
    }
    atomicParameter = new real_t[nzs];


    if (isNonlinear){
        nonlinearMat = new real_t*[np1+1]; // multiply matrix by f lnLc to get p*nu_s on p flux grid
        for (len_t i = 0; i<np1+1; i++){
            nonlinearMat[i] = new real_t[np1];
        }
        const real_t *p = mg->GetP1();
        trapzWeights = new real_t[np1];
        for (len_t i = 1; i<np1-1; i++){
            trapzWeights[i] = (p[i+1]-p[i-1])/2;
        }
    }
    
}


void CollisionFrequency::DeallocatePartialQuantities(){
    DeallocateGSL();
    if (nbound != nullptr){
        delete [] nbound;
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

    if(!buildOnlyF1F2){
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
    }
    if (screenedTerm_f1 != nullptr){
        delete [] screenedTerm_f1;
        delete [] screenedTerm_f2;
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
        for(len_t i = 0; i<np1+1;i++){
            delete [] nonlinearMat[i];
        }
        delete [] nonlinearMat;
        delete [] trapzWeights;
    }
}


/**
 * Initializes a GSL workspace for each radius (used for relativistic test particle operator evaluation),
 * using a T_cold-dependent fixed quadrature. 
 */
void CollisionFrequency::InitializeGSLWorkspace(){
 /** 
  * (consider using a single regular dynamic quadrature instead as the integral is somewhat tricky, 
  * since in the limit p/mc -> 0 the integral is sharply peaked at p_min -- goes as int 1/sqrt(x) dx,0,inf --
  * and may be challenging to resolve using a fixed point quadrature)
  */
    DeallocateGSL();
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    
    gsl_w = new gsl_integration_fixed_workspace*[nr];
    const real_t lowerLim = 0; // integrate from 0 to inf
    const gsl_integration_fixed_type *T = gsl_integration_fixed_laguerre;
    const len_t Npoints = 20; // play around with this number -- may require larger, or even sufficient with lower
    const real_t alpha = 0.0;
    real_t b;
    real_t Theta;
    for (len_t ir = 0; ir<nr; ir++){
        Theta = T_cold[ir]/Constants::mc2inEV;
        b = 1/Theta;
        gsl_w[ir] = gsl_integration_fixed_alloc(T, Npoints, lowerLim, b, alpha, 0.0);
    }
}


void CollisionFrequency::DeallocateGSL(){
    if (this->gsl_w == nullptr) // not sure if this works -- does freeing make it a nullptr?
        return;

    for (len_t ir=0; ir<this->nr; ir++)
        gsl_integration_fixed_free(gsl_w[ir]);
}


