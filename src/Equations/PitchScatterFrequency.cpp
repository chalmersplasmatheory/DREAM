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
                : CollisionFrequency(g,u,ih,lnLee,lnLei,mgtype,cqset){
}





/**
 * Evaluates the collision frequency at radial grid point ir and momentum p,
 * neglecting any contribution from the nonlinear collision operator
 */
real_t PitchScatterFrequency::evaluateAtP(len_t ir, real_t p){
    real_t *ncold = unknowns->GetUnknownData(id_ncold);
    real_t ntarget = ncold[ir];
    if (isNonScreened)
        ntarget += nbound[ir];

    real_t preFact = evaluatePreFactorAtP(p); 
    real_t lnLee = lnLambdaEE->evaluateAtP(ir,p);
    real_t lnLei = lnLambdaEI->evaluateAtP(ir,p);
    
    real_t collQty = lnLee * evaluateElectronTermAtP(ir,p) * ntarget;

    len_t ind;
    for(len_t iz = 0; iz<nZ; iz++){
        for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
            ind = ionIndex[iz][Z0];
            collQty += lnLei * evaluateIonTermAtP(iz,Z0,p) * ionDensities[ir][ind];
        }
    }
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

real_t PitchScatterFrequency::evaluateScreenedTermAtP(len_t iz, len_t Z0, real_t p){
    len_t ind = ionIndex[iz][Z0];
    len_t Z = Zs[iz];
    real_t a = atomicParameter[ind];
    real_t x = p*a*sqrt(p*a); 
    return 2.0/3.0 * ((Z*Z-Z0*Z0)*log(1+x) - (Z-Z0)*(Z-Z0)*x/(1+x) );
}

real_t PitchScatterFrequency::evaluateIonTermAtP(len_t iz, len_t Z0, real_t /*p*/){
    if (isNonScreened)
        return Zs[iz]*Zs[iz];
    else
        return Z0*Z0;
}

/**
 * Helper function to calculate a partial contribution to evaluateAtP
 */
real_t PitchScatterFrequency::evaluateElectronTermAtP(len_t ir, real_t p){
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
real_t PitchScatterFrequency::GetAtomicParameter(len_t iz, len_t Z0){
    len_t Z = ionHandler->GetZ(iz);
    // Fetch DFT-calculated value from table if it exists:
    for (len_t n=0; n<ionSizeAj_len; n++)
        if( Z==ionSizeAj_Zs[n] && (Z0==ionSizeAj_Z0s[n]) )
            return 2/Constants::alpha*ionSizeAj_data[n];

    // If DFT-data is missing, use Kirillov's model:
    return 2/Constants::alpha * pow(9*M_PI,1./3) / 4 * pow(Z-Z0,2./3) / Z;

}


void PitchScatterFrequency::GetNColdPartialContribution(len_t fluxGridMode, real_t *&partQty){
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
void PitchScatterFrequency::GetNiPartialContribution(len_t fluxGridMode, real_t *&partQty){
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

void PitchScatterFrequency::GetNonlinearPartialContribution(real_t *&partQty){
    GetNonlinearPartialContribution(lnLambdaEE->GetLnLambdaC(),partQty);
}

void PitchScatterFrequency::GetNiPartialContribution(real_t **nColdTerm, real_t *ionTerm, real_t *screenedTerm, real_t *preFactor, real_t *const* lnLee,  real_t *const* lnLei, len_t nr, len_t np1, len_t np2, real_t *&partQty){
    if(partQty==nullptr){
        partQty = new real_t[nzs*np1*np2*nr];
    }

    len_t pind, pindStore, indZ;
    for(len_t i = 0; i<np1; i++){
        for(len_t j = 0; j<np2; j++){
            pind = np1*j+i;
            if(isPXiGrid)
                pindStore = i;
            else
                pindStore = pind;
            
            for(len_t ir = 0; ir<nr; ir++){
                for(len_t iz=0; iz<nZ; iz++){
                    for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                        indZ = ionIndex[iz][Z0]; 
                        partQty[indZ*nr*np1*np2 + ir*np1*np2 + pind] = ionTerm[pindStore]*preFactor[pindStore]*lnLei[ir][pind];
                    }
                }
            }
        }
    }

    if(isNonScreened){
        real_t *electronContribution = new real_t[nr*np1*np2];
        GetNColdPartialContribution(nColdTerm,preFactor,lnLee,nr,np1,np2,electronContribution);
        for(len_t i = 0; i<np1; i++){
            for(len_t j = 0; j<np2; j++){
                pind = np1*j+i;
                if(isPXiGrid)
                    pindStore = i;
                else
                    pindStore = pind;
                
                for(len_t ir = 0; ir<nr; ir++){
                    for(len_t iz=0; iz<nZ; iz++){
                        for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                            indZ = ionIndex[iz][Z0]; 
                            partQty[indZ*nr*np1*np2 + ir*np1*np2 + pind] += (Zs[iz]-Z0)*electronContribution[ir*np1*np2 + pind];
                        }
                    }
                }
            }
        }
        delete [] electronContribution;
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

void PitchScatterFrequency::GetNColdPartialContribution(real_t **nColdTerm,real_t *preFactor, real_t *const* lnLee, len_t nr, len_t np1, len_t np2, real_t *&partQty){
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

void PitchScatterFrequency::GetNonlinearPartialContribution(const real_t* lnLc, real_t *&partQty){
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

