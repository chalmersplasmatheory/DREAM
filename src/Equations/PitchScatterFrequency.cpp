/**
 * Implementation of class which handles the pitch-angle scattering frequency nu_D.
 * Contains a pointer to an UnknownQuantityHandler; when plasma parameters have changed,
 * update frequency with nu_D->Rebuild();
 * When rebuilding grid, call nu_D->GridRebuilt().
 * evaluateAtP contains the cleanest implementation of the calculation of the frequency,
 * whereas AssembleQuantity is an optimised version which is used by the CollisionQuantityHandler.
 * 
 * TODO: np2_store should have +1 for _f2 terms when not pxigrid.
 */
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/FVMException.hpp"
#include <string>

using namespace DREAM;

/**
 * Effective 
 */
const len_t PitchScatterFrequency::ionSizeAj_len = 55; 
const real_t PitchScatterFrequency::ionSizeAj_data[ionSizeAj_len] = { 0.631757734322417, 0.449864664424796, 0.580073385681175, 0.417413282378673, 0.244965367639212, 0.213757911761448, 0.523908484242040, 0.432318176055981, 0.347483799585738, 0.256926098516580, 0.153148466772533, 0.140508604177553, 0.492749302776189, 0.419791849305259, 0.353418389488286, 0.288707775999513, 0.215438905215275, 0.129010899184783, 0.119987816515379, 0.403855887938967, 0.366602498048607, 0.329462647492495, 0.293062618368335, 0.259424839110224, 0.226161504309134, 0.190841656429844, 0.144834685411878, 0.087561370494245, 0.083302176729104, 0.351554934261205, 0.328774241757188, 0.305994557639981, 0.283122417984972, 0.260975850956140, 0.238925715853581, 0.216494264086975, 0.194295316086760, 0.171699132959493, 0.161221485564969, 0.150642403738712, 0.139526182041846, 0.128059339783537, 0.115255069413773, 0.099875435538094, 0.077085983503479, 0.047108093547224, 0.045962185039177, 0.235824746357894, 0.230045911002090, 0.224217341261303, 0.215062179624586, 0.118920957451653, 0.091511805821898, 0.067255603181663, 0.045824624741631 };
const real_t PitchScatterFrequency::ionSizeAj_Zs[ionSizeAj_len] = { 2, 2, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 54, 54, 54, 74, 74, 74, 74, 74 };
const real_t PitchScatterFrequency::ionSizeAj_Z0s[ionSizeAj_len] = { 0, 1, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 1, 2, 3, 0, 30, 40, 50, 60 };

/**
 * Constructor
 */
PitchScatterFrequency::PitchScatterFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                CoulombLogarithm *lnLei, CoulombLogarithm *lnLee,
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantityHandler::collqtyhand_settings *cqset)
                : CollisionQuantity(g,u,ih,mgtype,cqset){
    lnLambdaEI = lnLei;
    lnLambdaEE = lnLee;
}



/**
 * Evaluates the collision frequency at radial grid point ir and momentum p,
 * neglecting any contribution from the nonlinear collision operator
 */
real_t PitchScatterFrequency::evaluateAtP(len_t ir, real_t p){
    real_t gamma = sqrt(1+p*p);
    real_t *ncold = unknowns->GetUnknownData(id_ncold);
    real_t ntarget = ncold[ir];
    if (isNonScreened)
        ntarget += nbound[ir];

    real_t collQty;
    real_t preFact = evaluatePrefactorAtP(p,gamma)/gamma; 
    const len_t *Zs = ionHandler->GetZs();
    collQty = lnLambdaEE->evaluateAtP(ir,p) * evaluateGColdAtP(ir,p) * ntarget;
    for(len_t iz = 0; iz<nZ; iz++){
        for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
            collQty +=  lnLambdaEI->evaluateAtP(ir,p) * Z0*Z0 * ionHandler->GetIonDensity(ir,iz,Z0);
        }
    }
    if(isPartiallyScreened){
        for(len_t iz = 0; iz<nZ; iz++){
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                collQty +=  evaluateGKirillovAtP(iz,Z0,p) * ionHandler->GetIonDensity(ir,iz,Z0);
            }
        }
    }else if(isNonScreened){
        for(len_t iz = 0; iz<nZ; iz++){
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                collQty += (Zs[iz]*Zs[iz]-Z0*Z0)*lnLambdaEI->evaluateAtP(ir,p)*ionHandler->GetIonDensity(ir,iz,Z0);
            }
        }
    }
    collQty *= preFact;


    return collQty;
}
real_t PitchScatterFrequency::evaluateGKirillovAtP(len_t iz, len_t Z0, real_t p){
    len_t ind = ionHandler->GetIndex(iz,Z0);
    const len_t Z = ionHandler->GetZ(iz); 
    real_t a = ionEffectiveSize[ind];
    real_t x = p*a*sqrt(p*a); 
    return 2.0/3.0 * ((Z*Z-Z0*Z0)*log(1+x) - (Z-Z0)*(Z-Z0)*x/(1+x) );
                
}

/**
 * Helper function to calculate a partial contribution to evaluateAtP
 */
real_t PitchScatterFrequency::evaluateGColdAtP(len_t ir, real_t p){
    if (collQtySettings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL){
        if(p==0)
            return ReallyLargeNumber;
        else 
            return 1;
    } else if (collQtySettings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
        if(p==0)
            return 0;
        real_t p2 = p*p;
        real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
        real_t Theta,M;
        real_t gamma = sqrt(1+p2);
        Theta = T_cold[ir] / Constants::mc2inEV;
        M = 0;
        M += (p2*gamma*gamma + Theta*Theta)*evaluatePsi0(ir,p);
        M += Theta*(2*p2*p2 - 1)*evaluatePsi1(ir,p);
        M += gamma*Theta * ( 1 + Theta*(2*p2-1)*p*exp( -(gamma-1)/Theta ) );
        M /= gamma*gamma*p*p*evaluateExp1OverThetaK(Theta,2.0);
        return  M;
    } else {
        throw FVM::FVMException("Invalid collfreq_mode.");
        return -1;
    }
}

/**
 * Rebuilds the partial contributions to the collision frequency that depend on density and temperature
 */
void PitchScatterFrequency::RebuildPlasmaDependentTerms(){
    nbound = ionHandler->evaluateBoundElectronDensityFromQuasiNeutrality(nbound);
    if (!buildOnlyF1F2){
        setGCold(gCold,mg->GetP(),nr,np1,np2_store);
        setGCold(gCold_fr,mg->GetP(),nr+1,np1,np2_store);
    }
    setGCold(gCold_f1,mg->GetP_f1(),nr,np1+1,np2_store);
    setGCold(gCold_f2,mg->GetP_f2(),nr,np1,np2_store);
}

/**
 * Rebuilds partial contributions that only depend on the grid. If P-Xi grid, only store momentum
 * dependent quantities on size np1 array. 
 */
void PitchScatterFrequency::RebuildConstantTerms(){
    
    const len_t *Zs = ionHandler->GetZs();


    len_t ind;
    for(len_t iZ = 0; iZ<nZ; iZ++){
        for(len_t Z0=0; Z0<=Zs[iZ]; Z0++){
            ind = ionHandler->GetIndex(iZ,Z0);
            ionEffectiveSize[ind] = GetIonEffectiveSizeAj(iZ,Z0);
        }
    }

    if (!buildOnlyF1F2){
        setPreFactor(preFactor,mg->GetP(),mg->GetGamma(),np1,np2_store);
        setPreFactor(preFactor_fr,mg->GetP(),mg->GetGamma(),np1,np2_store);
        if(isPartiallyScreened){
            setGKirillov(gKirillov,mg->GetP(),np1,np2_store);
            setGKirillov(gKirillov_fr,mg->GetP(),np1,np2_store);
        }
    }
    setPreFactor(preFactor_f1,mg->GetP_f1(),mg->GetGamma_f1(),np1+1,np2_store);
    setPreFactor(preFactor_f2,mg->GetP_f2(),mg->GetGamma_f2(),np1,np2_store);
    if(isPartiallyScreened){
        setGKirillov(gKirillov_f1,mg->GetP_f1(),np1+1,np2_store);
        setGKirillov(gKirillov_f2,mg->GetP_f2(),np1,np2_store);
    }
    if(isNonlinear)
        calculateIsotropicNonlinearOperatorMatrix();

}

/**
 * Calculates and stores the partial contribution to the collision frequency from free electrons
 */
void PitchScatterFrequency::setGCold(real_t **&gCold, const real_t *pIn, len_t nr, len_t np1, len_t np2){
    real_t p;
    len_t pind;
    // Depending on setting, set nu_s to superthermal or full formula (with maxwellian)
    for(len_t i=0;i<np1;i++){
        for(len_t j=0;j<np2;j++){
            pind = np1*j+i;
            p = pIn[pind];
            for(len_t ir=0;ir<nr;ir++)
                gCold[ir][pind] = evaluateGColdAtP(ir,p);
        }
    }
    

}

/**
 * Calculates and stores the partial contribution to the collision frequency from bound electrons using the Bethe formula.
 */
void PitchScatterFrequency::setGKirillov(real_t *&gKirillov, const real_t *pIn, len_t np1, len_t np2){
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
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                    ind = ionHandler->GetIndex(iz,Z0);
                    gKirillov[ind*np1*np2 + pind] = evaluateGKirillovAtP(iz,Z0,p);
                }
            }
        }
    }
}

/**
 * Calculates and stores the momentum-dependent prefactor to the collision frequency.
 */
void PitchScatterFrequency::setPreFactor(real_t *&preFactor, const real_t *pIn, const real_t *gammaIn, len_t np1, len_t np2){
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
                PF = evaluatePrefactorAtP(p,gamma)/gamma;
            preFactor[ind] = PF;
        }
    }
}




// Calculates Rosenbluth potential matrices defined such that when they are muliplied
// by the f_hot distribution vector, yields the three collision frequencies.
// XXX assuming for now same grid at all radii, and hot tail grid (PXi, nxi=1). 
void PitchScatterFrequency::calculateIsotropicNonlinearOperatorMatrix(){
    


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
        p2 = p[0]*p[0];
        nonlinearMat[i][0] = (4*M_PI/3) * constPreFactor / p_f[i]*( (p[1]-p[0])/2*(3-p2/p2f) + p[0]*(1-p2/(5*p2f) ))*p2/p2f;
        for (len_t ip = 1; ip < i-1; ip++){
            p2 = p[ip]*p[ip];
            nonlinearMat[i][ip] = (4*M_PI/3) * constPreFactor / p_f[i] * trapzWeights[ip]*p2/p2f *(3-p2/p2f) ;
        } 
        p2 = p[i-1]*p[i-1];
        weightsIm1 = (p[i-1]-p[i-2])/2 + (p_f[i]-p[i-1])/(p[i]-p[i-1])*( (2*p[i]-p_f[i]-p[i-1])/2 );
        nonlinearMat[i][i-1] = (4*M_PI/3) * constPreFactor / p_f[i] * weightsIm1*p2/p2f *(3-p2/p2f) ;
        p2 = p[i]*p[i];
        weightsI = (p_f[i]-p[i-1])*(p_f[i]-p[i-1])/(p[i]-p[i-1]);
        nonlinearMat[i][i] = (4*M_PI/3) * constPreFactor / p_f[i] * weightsI*p2/p2f *(3-p2/p2f) ;

        // add contribution from p'>p terms near p'=p
        p2 = p[i-1]*p[i-1];
        weightsIm1 = (1.0/2)*(p[i]-p_f[i])*(p[i]-p_f[i])/(p[i]-p[i-1]);
        nonlinearMat[i][i-1] += (8*M_PI/3) * constPreFactor / p_f[i] * weightsIm1*p[i-1]/p2f;
        p2 = p[i]*p[i];
        weightsI = (p[i+1]-p[i])/2 + (1.0/2)*(p[i]-p_f[i])*(p_f[i]+p[i]-2*p[i-1])/(p[i]-p[i-1]);
        nonlinearMat[i][i] += (8*M_PI/3) * constPreFactor * weightsI * p[i]/p2f;


        for (len_t ip = i+1; ip < np1-1; ip++){
            nonlinearMat[i][ip] = (8*M_PI/3) * constPreFactor * trapzWeights[ip]*p[ip]/p2f ;
        } 
        real_t weightsEnd = (p[np1-1]-p[np1-2])/2;
        nonlinearMat[i][np1-1] = (8*M_PI/3) * constPreFactor * weightsEnd*p[np1-1]/p2f ;
        

    }
}

/**
 * Returns the mean excitation energy of ion species with index iz and charge state Z0.
 * Currently only for ions we have actual data for -- should investigate approximate models for other species
 */
real_t PitchScatterFrequency::GetIonEffectiveSizeAj(len_t iz, len_t Z0){
    len_t Z = ionHandler->GetZ(iz);
    // Fetch DFT-calculated value from table if it exists:
    for (len_t n=0; n<ionSizeAj_len; n++)
        if( Z==ionSizeAj_Zs[n] && (Z0==ionSizeAj_Z0s[n]) )
            return 2/Constants::alpha*ionSizeAj_data[n];

    // If DFT-data is missing, use Kirillov's model:
    return 2/Constants::alpha * pow(9*M_PI,1./3) / 4 * pow(Z-Z0,2./3) / Z;

}


/**
 * The following methods returns the partial contribution to the collision frequency from the unknownQuantity
 * with index id_unknown (corresponding either to "n_cold", "n_i" or "f_hot" (for nonlinear operator)). 
 * I.e., returns the partial derivative of nu_s with respect to the corresponding density (neglecting logarithmic
 * variations in lnLambda), for the evaluation of the Jacobian in the Newton method.  
 */
void PitchScatterFrequency::GetPartialContribution(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty){
    len_t pind;
    if(isPXiGrid)
        pind = i;
    else 
        pind = np1*j+i;
    evaluatePartialContribution(id_unknown,preFactor,gKirillov,gCold,lnLambdaEI->GetValue(ir,i,j),lnLambdaEE->GetValue(ir,i,j),ir,pind,np1,partQty);
}
void PitchScatterFrequency::GetPartialContribution_fr(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty){
    len_t pind;
    if(isPXiGrid)
        pind = i;
    else 
        pind = np1*j+i;
    evaluatePartialContribution(id_unknown,preFactor_fr,gKirillov_fr,gCold_fr,lnLambdaEI->GetValue_fr(ir,i,j),lnLambdaEE->GetValue_fr(ir,i,j),ir,pind,np1,partQty);
}
void PitchScatterFrequency::GetPartialContribution_f1(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty){
    len_t pind;
    if(isPXiGrid)
        pind = i;
    else 
        pind = (np1+1)*j+i;
    evaluatePartialContribution(id_unknown,preFactor_f1,gKirillov_f1,gCold_f1,lnLambdaEI->GetValue_f1(ir,i,j),lnLambdaEE->GetValue_f1(ir,i,j),ir,pind,np1+1,partQty);
}
void PitchScatterFrequency::GetPartialContribution_f2(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty){
    len_t pind;
    if(isPXiGrid)
        pind = i;
    else 
        pind = np1*j+i;
    evaluatePartialContribution(id_unknown,preFactor_f2,gKirillov_f2,gCold_f2,lnLambdaEI->GetValue_f2(ir,i,j),lnLambdaEE->GetValue_f2(ir,i,j),ir,pind,np1,partQty);
}

void PitchScatterFrequency::evaluatePartialContribution(len_t id_unknown, real_t *preFactor, real_t *gKirillov, real_t **gCold, real_t lnLei, real_t lnLee, len_t ir, len_t pind,len_t np1,real_t *&partQty){
    
    if (gKirillov == nullptr)
        FVM::FVMException("Quantities have not been built on the distribution grid. Set buildOnlyF1F2=false.");

    if(id_unknown == id_ni)
        GetPartialContributionNi(preFactor[pind],gKirillov,gCold[ir][pind],lnLei,lnLee,pind,np1, partQty);
    else if (id_unknown == id_ncold){
        GetPartialContributionNCold(gCold[ir][pind],preFactor[pind],lnLee, partQty);
    } else if(id_unknown == id_fhot)
        GetPartialContributionNonlinear(lnLambdaEE->GetLnLambdaC(ir),pind,np1, partQty);
}

void PitchScatterFrequency::GetPartialContributionNCold(real_t gCold, real_t preFactor, real_t lnLee, real_t *&partQty){
    if(partQty==nullptr)
        partQty = new real_t[1];
    partQty[0] = gCold * preFactor * lnLee;
}

void PitchScatterFrequency::GetPartialContributionNi(real_t preFactor, real_t *gKirillov, real_t gCold, const real_t lnLei, const real_t lnLee, len_t pind, len_t np1,real_t *&partQty){
    if(partQty==nullptr)
        partQty = new real_t[nzs];

    len_t indZ;
    const len_t *Zs = ionHandler->GetZs();

    for(len_t iz = 0; iz<nZ; iz++){
        for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
            indZ = ionHandler->GetIndex(iz,Z0);
            partQty[indZ] = Z0*Z0*preFactor*lnLei ;
        }
    }
    if(collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_COMPLETELY_SCREENED){
        return;
    } else if (collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED){

        real_t *boundElectronContribution = new real_t[1];
        GetPartialContributionNCold(gCold, preFactor, lnLee, boundElectronContribution);
        for(len_t iz = 0; iz<nZ; iz++){
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                indZ = ionHandler->GetIndex(iz,Z0);
                partQty[indZ] +=  (Zs[iz]-Z0)*boundElectronContribution[0] + preFactor * lnLei * ( Zs[iz]*Zs[iz]-Z0*Z0 );
            }
        }
        delete [] boundElectronContribution;
    } else if (collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED){
        for(len_t iz = 0; iz<nZ; iz++){
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                indZ = ionHandler->GetIndex(iz,Z0);
                partQty[indZ] +=  preFactor * gKirillov[indZ*np1*np2_store + pind];
            }
        }
    }
}

void PitchScatterFrequency::GetPartialContributionNonlinear(real_t lnLc, len_t pind, len_t np1,real_t *&partQty){
    if(partQty==nullptr)
        partQty = new real_t[nr*np1];
    for(len_t ir=0;ir<nr;ir++)
        for(len_t ip=0; ip<np1; ip++)
            partQty[np1*ir+ip] = lnLc*nonlinearMat[pind][ip];
            
}



/**
 * Allocates quantities which will be used in the calculation of the collision frequency.
 */
void PitchScatterFrequency::AllocatePartialQuantities(){
    
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
            gKirillov    = new real_t[nzs*np1*np2_store];
            gKirillov_fr = new real_t[nzs*np1*np2_store];
        }
        AllocateGColdFunc(gCold,nr,np1,np2_store);
        AllocateGColdFunc(gCold_fr,nr+1,np1,np2_store);
    }
    if(isPartiallyScreened){
        gKirillov_f1 = new real_t[nzs*(np1+1)*np2_store];
        gKirillov_f2 = new real_t[nzs*np1*np2_store];
    }
    AllocateGColdFunc(gCold_f1,nr,np1+1,np2_store);
    AllocateGColdFunc(gCold_f2,nr,np1,np2_store);
    ionEffectiveSize = new real_t[nzs];
    
}

void PitchScatterFrequency::DeallocatePartialQuantities(){
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
        DeallocateGColdFunc(gCold, nr);
        DeallocateGColdFunc(gCold_fr, nr+1);
        if (gKirillov != nullptr)
            delete [] gKirillov;
        if (gKirillov_fr != nullptr)
            delete [] gKirillov_fr;
    }
    if (gKirillov_f1 != nullptr)
        delete [] gKirillov_f1;
    if (gKirillov_f2 != nullptr)
        delete [] gKirillov_f2;
    DeallocateGColdFunc(gCold_f1, nr);
    DeallocateGColdFunc(gCold_f2, nr);

    if(ionEffectiveSize!=nullptr)
        delete [] ionEffectiveSize;

}


void PitchScatterFrequency::AllocateGColdFunc(real_t **&gCold,len_t nr,len_t np1, len_t np2){
    gCold = new real_t*[nr];
    for(len_t ir=0;ir<nr;ir++)
        gCold[ir] = new real_t[np1*np2];
}

void PitchScatterFrequency::DeallocateGColdFunc(real_t **&gCold, len_t nr){
    if(gCold != nullptr){
        for(len_t ir; ir<nr; ir++)
            delete [] gCold[ir];
         delete [] gCold;
    }
}

