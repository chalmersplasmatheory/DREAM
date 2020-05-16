/**
 * Implementation of class which handles the slowing down frequency nu_s.
 * Contains a pointer to an UnknownQuantityHandler; when plasma parameters have changed,
 * update frequency with nu_s->Rebuild();
 * When rebuilding grid, call nu_s->GridRebuilt().
 * Can return nu_s in two different ways that should be identical (up to some machine precision error)
 * 1:
 * p = mg->GetP(i,j)
 * val=nuS->evaluateAtP(ir,p)
 * 2:
 * val=nuS->GetValue(ir,i,j)
 *
 * evaluateAtP contains the cleanest implementation of the calculation of the frequency,
 * whereas AssembleQuantity is an optimised version which is used by the CollisionQuantityHandler.
 * 
 * TODO: we should design a test that verifies that these two are equal for each 
 * collfreq_type, collfreq_mode setting and for each collision frequency, and 
 * also benchmarks to a separate implementation (say, the CODE one).
 * Fix np2_store and _f2 terms for non-pxi grids
 */
#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/FVMException.hpp"
#include <string>

using namespace DREAM;

/**
 * Mean excitation energy atomic data for ions from 
 *   Sauer, S.P., Oddershede, J. and Sabin, J.R., 2015. 
 *   The mean excitation energy of atomic ions. 
 *   In Advances in Quantum Chemistry (Vol. 71, pp. 29-40). 
 *   Academic Press.
 */
const len_t SlowingDownFrequency::meanExcI_len = 40;
const real_t SlowingDownFrequency::meanExcI_data[meanExcI_len] = {2.9295e-5, 8.3523e-05, 1.1718e-04, 6.4775e-05, 2.1155e-04, 2.6243e-04, 1.2896e-04, 1.8121e-04, 
        2.6380e-04, 4.1918e-04, 9.5147e-04, 0.0011, 2.6849e-04, 3.2329e-04, 3.8532e-04, 4.6027e-04, 5.5342e-04, 
        6.9002e-04, 9.2955e-04, 0.0014, 0.0028, 0.0029, 3.6888e-04, 4.2935e-04, 4.9667e-04, 5.7417e-04, 6.6360e-04, 
        7.7202e-04, 9.0685e-04, 0.0011, 0.0014, 0.0016, 0.0017, 0.0019, 0.0022, 0.0027, 0.0035, 0.0049, 0.0092, 0.0095};
const real_t SlowingDownFrequency::meanExcI_Zs[meanExcI_len] = {1, 2, 2, 3, 3, 3, 6, 6, 6, 6, 6, 6, 10, 10, 10, 10, 10, 10, 
        10, 10, 10, 10, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18};
const real_t SlowingDownFrequency::meanExcI_Z0s[meanExcI_len] = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};


/**
 * Constructor
 */
SlowingDownFrequency::SlowingDownFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                CoulombLogarithm *lnLee,
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantityHandler::collqtyhand_settings *cqset)
                : CollisionQuantity(g,u,ih,mgtype,cqset){
    lnLambdaEE = lnLee;
}


/**
 * Evaluates the collision frequency at radial grid point ir and momentum p,
 * neglecting any contribution from the nonlinear collision operator
 */
real_t SlowingDownFrequency::evaluateAtP(len_t ir, real_t p){
    real_t gamma = sqrt(1+p*p);
    real_t *ncold = unknowns->GetUnknownData(id_ncold);
    real_t ntarget = ncold[ir];
    if (isNonScreened)
        ntarget += nbound[ir];

    real_t collQty;
    real_t preFact = evaluatePrefactorAtP(p,gamma); 
    const len_t *Zs = ionHandler->GetZs();
    collQty = lnLambdaEE->evaluateAtP(ir,p) * evaluateHColdAtP(ir,p) * ntarget;
    if(isPartiallyScreened){
        
        for(len_t iz = 0; iz<nZ; iz++){
            for(len_t Z0=0; Z0<Zs[iz]; Z0++){
                collQty +=  evaluateHBetheAtP(iz,Z0,p) * ionHandler->GetIonDensity(ir,iz,Z0);
            }
        }
    }
    collQty *= preFact;


    return collQty;
}


real_t SlowingDownFrequency::evaluateHBetheAtP(len_t iz, len_t Z0, real_t p){
    len_t Z = ionHandler->GetZ(iz); 
    len_t ind = ionHandler->GetIndex(iz,Z0);
    real_t gamma = sqrt(1+p*p);
    real_t beta = p/gamma;
    real_t h = p*sqrt(gamma-1)/meanExcitationEnergy[ind];
    real_t nBound = Z - Z0;
    return nBound*(log(1+pow(h,kInterpolate))/kInterpolate - beta*beta) ;

}

/**
 * Helper function to calculate a partial contribution to evaluateAtP
 */
real_t SlowingDownFrequency::evaluateHColdAtP(len_t ir, real_t p){
    if (collQtySettings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL){
        if(p==0)
            return 1e50;
        else 
            return 1;
    } else if (collQtySettings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
        if(p==0)
            return 0;
        real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
        real_t Theta,M;
        real_t gamma = sqrt(1+p*p);
        Theta = T_cold[ir] / Constants::mc2inEV;
        M = 0;
        M += gamma*gamma* evaluatePsi1(Theta,p) - Theta * evaluatePsi0(Theta,p);
        M +=  (Theta*gamma - 1) * p * exp( -(gamma-1)/Theta );
        M /= evaluateExp1OverThetaK(Theta,2.0);
        return  M  / (gamma*gamma);
    } else {
        throw FVM::FVMException("Invalid collfreq_mode.");
        return -1;
    }
}

/**
 * Rebuilds the partial contributions to the collision frequency that depend on density and temperature
 */
void SlowingDownFrequency::RebuildPlasmaDependentTerms(){
    nbound = ionHandler->evaluateBoundElectronDensityFromQuasiNeutrality(nbound);
    if (!buildOnlyF1F2){
        setHCold(hCold,mg->GetP(),mg->GetGamma(),nr,np1,np2_store);
        setHCold(hCold_fr,mg->GetP(),mg->GetGamma(),nr+1,np1,np2_store);
    }
    setHCold(hCold_f1,mg->GetP_f1(),mg->GetGamma_f1(),nr,np1+1,np2_store);
    setHCold(hCold_f2,mg->GetP_f2(),mg->GetGamma_f2(),nr,np1,np2_store);
}

/**
 * Rebuilds partial contributions that only depend on the grid. If P-Xi grid, only store momentum
 * dependent quantities on size np1 array. 
 */
void SlowingDownFrequency::RebuildConstantTerms(){
    const len_t *Zs = ionHandler->GetZs();

    len_t ind;
    for(len_t iZ = 0; iZ<nZ; iZ++){
        for(len_t Z0=0; Z0<=Zs[iZ]-1; Z0++){
            ind = ionHandler->GetIndex(iZ,Z0);
            meanExcitationEnergy[ind] = GetMeanExcitationEnergy(iZ,Z0);
        }
    }

    if (!buildOnlyF1F2){
        setPreFactor(preFactor,mg->GetP(),mg->GetGamma(),np1,np2_store);
        setPreFactor(preFactor_fr,mg->GetP(),mg->GetGamma(),np1,np2_store);
        if(isPartiallyScreened){
            setHBethe(hiBethe,mg->GetP(),np1,np2_store);
            setHBethe(hiBethe_fr,mg->GetP(),np1,np2_store);
        }
    }
    setPreFactor(preFactor_f1,mg->GetP_f1(),mg->GetGamma_f1(),np1+1,np2_store);
    setPreFactor(preFactor_f2,mg->GetP_f2(),mg->GetGamma_f2(),np1,np2_store);
    if(isPartiallyScreened){
        setHBethe(hiBethe_f1,mg->GetP_f1(),np1+1,np2_store);
        setHBethe(hiBethe_f2,mg->GetP_f2(),np1,np2_store);
    }
    if(isNonlinear)
        calculateIsotropicNonlinearOperatorMatrix();

}

/**
 * Calculates and stores the partial contribution to the collision frequency from free electrons
 */
void SlowingDownFrequency::setHCold(real_t **&hCold, const real_t *pIn, const real_t *gammaIn, len_t nr, len_t np1, len_t np2){
    real_t h,p;
    len_t pind;
    // Depending on setting, set nu_s to superthermal or full formula (with maxwellian)
    if (collQtySettings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL) {
        for(len_t i=0;i<np1;i++){
            for(len_t j=0;j<np2;j++){
                pind = np1*j+i;
                p = pIn[pind];
                if(p==0)
                    h = 1e50;
                else
                    h = 1;
                for(len_t ir=0; ir<nr; ir++)
                    hCold[ir][pind] = h;
            }
        }
    } else if (collQtySettings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
        real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
        real_t gamma,Theta,M;
        for(len_t i=0;i<np1;i++){
            for(len_t j=0;j<np2;j++){
                pind = np1*j+i;
                p = pIn[pind];
                gamma = gammaIn[pind];
                if(p==0)
                    for(len_t ir=0; ir<nr; ir++)
                        hCold[ir][np1*j+i] = 0;
                else {
                    for(len_t ir=0; ir<nr; ir++){
                        Theta = T_cold[ir] / Constants::mc2inEV;
                        M = 0;
                        M += gamma*gamma* evaluatePsi1(Theta,p) - Theta * evaluatePsi0(Theta,p);
                        M +=  (Theta*gamma - 1) * p * exp( -(gamma-1)/Theta );
                        M /= evaluateExp1OverThetaK(Theta,2.0);
                        hCold[ir][np1*j+i]  =  M  / (gamma*gamma);
                    }
                }
            }
        }
    } else
        throw NotImplementedException("Chosen collfreq_mode setting not yet supported.");

}



/**
 * Calculates and stores the partial contribution to the collision frequency from bound electrons using the Bethe formula.
 */
void SlowingDownFrequency::setHBethe(real_t *&hiBethe, const real_t *pIn, len_t np1, len_t np2){
    if (!(collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED))
        return;

    const len_t *Zs = ionHandler->GetZs();
    real_t p;
    len_t ind, pind;
    for(len_t i = 0; i<np1; i++){
        for (len_t j = 0; j<np2; j++){
            pind = np1*j+i;
            p = pIn[pind];
            for(len_t iz = 0; iz<nZ; iz++){
                for(len_t Z0=0; Z0<Zs[iz]; Z0++){
                    ind = ionHandler->GetIndex(iz,Z0);
                    hiBethe[ind*np1*np2 + pind] = evaluateHBetheAtP(iz,Z0,p);
                }
            }
        }
    }
}

/**
 * Calculates and stores the momentum-dependent prefactor to the collision frequency.
 */
void SlowingDownFrequency::setPreFactor(real_t *&preFactor, const real_t *pIn, const real_t *gammaIn, len_t np1, len_t np2){
    real_t p, gamma, PF;
    len_t ind;
    for(len_t i = 0; i<np1; i++){
        for (len_t j = 0; j<np2; j++){
            ind = np1*j+i;
            p = pIn[ind];
            gamma = gammaIn[ind];
            if(p==0)
                PF = ReallyLargeNumber; 
            else
                PF = evaluatePrefactorAtP(p,gamma);
            preFactor[ind] = PF;
        }
    }
}




// Calculates Rosenbluth potential matrices defined such that when they are muliplied
// by the f_hot distribution vector, yields the three collision frequencies.
// XXX assuming for now same grid at all radii, and hot tail grid (PXi, nxi=1). 
void SlowingDownFrequency::calculateIsotropicNonlinearOperatorMatrix(){

    if( !(isPXiGrid && (mg->GetNp2() == 1)) )
        throw NotImplementedException("Nonlinear collisions only implemented for hot tails (np2=1) and p-xi grid");

    
    const real_t *p_f = mg->GetP1_f();
    const real_t *p = mg->GetP1();


    // See doc/notes/theory.pdf appendix B for details on discretization of integrals;
    // uses a trapezoidal rule
    real_t p2, p2f;
    real_t weightsIm1, weightsI;
    for (len_t i = 1; i<np1+1; i++){
        p2f = p_f[i]*p_f[i];
        p2  = p[0]*p[0];
        nonlinearMat[i][0] = 4*M_PI/p_f[i] * constPreFactor*( (p[1]-p[0])/2 + p[0]/3 )*p2/p2f;
        for (len_t ip = 1; ip < i-1; ip++){
            p2 = p[ip]*p[ip];
            nonlinearMat[i][ip] = 4*M_PI/p_f[i] * constPreFactor* trapzWeights[ip]*p2/p2f;
        } 
        p2 = p[i-1]*p[i-1];
        weightsIm1 = (p[i-1]-p[i-2])/2 + (p_f[i]-p[i-1])/(p[i]-p[i-1])*( (2*p[i]-p_f[i]-p[i-1])/2 );
        nonlinearMat[i][i-1] = 4*M_PI/p_f[i] * constPreFactor * weightsIm1*p2/p2f ;
        p2 = p[i]*p[i];
        weightsI = (p_f[i]-p[i-1])*(p_f[i]-p[i-1])/(p[i]-p[i-1]);
        nonlinearMat[i][i]   = 4*M_PI/p_f[i] * constPreFactor * (1.0/2)* weightsI *p2/p2f;
    }
}


/**
 * Returns the mean excitation energy of ion species with index iz and charge state Z0.
 * Currently only for ions we have actual data for -- should investigate approximate models for other species
 */
real_t SlowingDownFrequency::GetMeanExcitationEnergy(len_t iz, len_t Z0){
    len_t Z = ionHandler->GetZ(iz);
    if(Z0==(Z-1)) // For hydrogenic ions the mean excitation energy is 14.9916*Z^2 eV 
        return Z*Z*14.9916 / Constants::mc2inEV;

    // Fetch value from table if it exists:
    for (len_t n=0; n<meanExcI_len; n++)
        if( (Z==meanExcI_Zs[n]) && (Z0==meanExcI_Z0s[n]) )
            return meanExcI_data[n];

    throw FVM::FVMException("Mean excitation energy for ion species: '%s' in charge state Z0 = " LEN_T_PRINTF_FMT " is missing.", ionHandler->GetName(iz).c_str(), Z0); 
}


/**
 * The following methods returns the partial contribution to the collision frequency from the unknownQuantity
 * with index id_unknown (corresponding either to "n_cold", "n_i" or "f_hot" (for nonlinear operator)). 
 * I.e., returns the partial derivative of nu_s with respect to the corresponding density (neglecting logarithmic
 * variations in lnLambda), for the evaluation of the Jacobian in the Newton method.  
 */
void SlowingDownFrequency::GetPartialContribution(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty){
    len_t pind;
    if(isPXiGrid)
        pind = i;
    else 
        pind = np1*j+i;
    evaluatePartialContribution(id_unknown,preFactor,hiBethe,hCold,lnLambdaEE->GetValue(ir,i,j),ir,pind,np1,partQty);
}
void SlowingDownFrequency::GetPartialContribution_fr(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty){
    len_t pind;
    if(isPXiGrid)
        pind = i;
    else 
        pind = np1*j+i;
    evaluatePartialContribution(id_unknown,preFactor_fr,hiBethe_fr,hCold_fr,lnLambdaEE->GetValue_fr(ir,i,j),ir,pind,np1,partQty);
}
void SlowingDownFrequency::GetPartialContribution_f1(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty){
    len_t pind;
    if(isPXiGrid)
        pind = i;
    else 
        pind = (np1+1)*j+i;
    evaluatePartialContribution(id_unknown,preFactor_f1,hiBethe_f1,hCold_f1,lnLambdaEE->GetValue_f1(ir,i,j),ir,pind,np1+1,partQty);
}
void SlowingDownFrequency::GetPartialContribution_f2(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty){
    len_t pind;
    if(isPXiGrid)
        pind = i;
    else 
        pind = np1*j+i;
    evaluatePartialContribution(id_unknown,preFactor_f2,hiBethe_f2,hCold_f2,lnLambdaEE->GetValue_f2(ir,i,j),ir,pind,np1,partQty);
}

void SlowingDownFrequency::evaluatePartialContribution(len_t id_unknown, real_t *preFactor, real_t *hiBethe, real_t **hCold, real_t lnLee, len_t ir, len_t pind,len_t np1, real_t *&partQty){
    
    if (hiBethe == nullptr)
        FVM::FVMException("Quantities have not been built on the distribution grid. Set buildOnlyF1F2=false.");

    if(id_unknown == id_ni)
        GetPartialContributionNi(preFactor[pind],hiBethe,hCold[ir][pind],lnLee,pind,np1, partQty);
    else if (id_unknown == id_ncold){
        if(partQty==nullptr)
            partQty = new real_t[1];
        partQty[0] = hCold[ir][pind] * preFactor[pind] * lnLee;
    } else if(id_unknown == id_fhot)
        GetPartialContributionNonlinear(lnLambdaEE->GetLnLambdaC(ir),pind, np1,partQty);
}

void SlowingDownFrequency::GetPartialContributionNi(real_t preFactor, real_t *hiBethe, real_t hCold, const real_t lnLee, len_t pind, len_t np1, real_t *&partQty){
    if(partQty==nullptr)
        partQty = new real_t[nzs];
    for(len_t it = 0; it<nzs; it++){
        partQty[it] = 0;
    }
    len_t indZ;
    const len_t *Zs = ionHandler->GetZs();

    if(collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_COMPLETELY_SCREENED){
        return;
    } else if (collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED){
        for(len_t iz = 0; iz<nZ; iz++){
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                indZ = ionHandler->GetIndex(iz,Z0);
                partQty[indZ] +=  preFactor * hCold * lnLee * (Zs[iz]-Z0);
            }
        }
    } else if (collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED){
        for(len_t iz = 0; iz<nZ; iz++){
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                indZ = ionHandler->GetIndex(iz,Z0);
                partQty[indZ] +=  preFactor * hiBethe[indZ*np1*np2_store + pind];
            }
        }
    }
}

void SlowingDownFrequency::GetPartialContributionNonlinear(real_t lnLc, len_t pind, len_t np1, real_t *&partQty){
    if(partQty==nullptr)
        partQty = new real_t[nr*np1];
    for(len_t ir=0;ir<nr;ir++)
        for(len_t ip=0; ip<np1; ip++)
            partQty[np1*ir+ip] = lnLc*nonlinearMat[pind][ip];
            
}



/**
 * Allocates quantities which will be used in the calculation of the collision frequency.
 */
void SlowingDownFrequency::AllocatePartialQuantities(){
   
    DeallocatePartialQuantities();
    nbound = new real_t[nr];
    if(!buildOnlyF1F2){
        preFactor    = new real_t[np1*np2_store];
        preFactor_fr = new real_t[np1*np2_store];
    }
    preFactor_f1 = new real_t[(np1+1)*np2_store];
    preFactor_f2 = new real_t[np1*np2_store];
    
    if(!buildOnlyF1F2){
        if(isPartiallyScreened){
            hiBethe    = new real_t[nzs*np1*np2_store];
            hiBethe_fr = new real_t[nzs*np1*np2_store];
        }
        AllocateHColdFunc(hCold,nr,np1,np2_store);
        AllocateHColdFunc(hCold_fr,nr+1,np1,np2_store);
    }
    if(isPartiallyScreened){
        hiBethe_f1 = new real_t[nzs*(np1+1)*np2_store];
        hiBethe_f2 = new real_t[nzs*np1*np2_store];
    }
    AllocateHColdFunc(hCold_f1,nr,np1+1,np2_store);
    AllocateHColdFunc(hCold_f2,nr,np1,np2_store);
    meanExcitationEnergy = new real_t[nzs];

    
}

void SlowingDownFrequency::DeallocatePartialQuantities(){
    if (nbound != nullptr)
        delete [] nbound;
    if(preFactor!=nullptr){
        delete [] preFactor;
        delete [] preFactor_fr;
    }   
    if(preFactor_f1 != nullptr){
        delete [] preFactor_f1;
        delete [] preFactor_f2;
    }

    if(!buildOnlyF1F2){
        DeallocateHColdFunc(hCold, nr);
        DeallocateHColdFunc(hCold_fr, nr+1);
        if (hiBethe != nullptr)
            delete [] hiBethe;
        if (hiBethe_fr != nullptr)
            delete [] hiBethe_fr;
    }
    if (hiBethe_f1 != nullptr)
        delete [] hiBethe_f1;
    if (hiBethe_f2 != nullptr)
        delete [] hiBethe_f2;
    DeallocateHColdFunc(hCold_f1, nr);
    DeallocateHColdFunc(hCold_f2, nr);

    if(meanExcitationEnergy != nullptr)
        delete [] meanExcitationEnergy;

}


void SlowingDownFrequency::AllocateHColdFunc(real_t **&hCold,len_t nr,len_t np1, len_t np2){
    hCold = new real_t*[nr];
    for(len_t ir=0;ir<nr;ir++)
        hCold[ir] = new real_t[np1*np2];
}

void SlowingDownFrequency::DeallocateHColdFunc(real_t **&hCold, len_t nr){
    if(hCold != nullptr){
        for(len_t ir; ir<nr; ir++)
            delete [] hCold[ir];
         delete [] hCold;
    }
}

