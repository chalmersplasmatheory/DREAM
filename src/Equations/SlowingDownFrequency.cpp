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
                CoulombLogarithm *lnLee,CoulombLogarithm *lnLei,
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantityHandler::collqtyhand_settings *cqset)
                : CollisionFrequency(g,u,ih,lnLee,lnLei,mgtype,cqset){
    hasIonTerm = false;
}

/**
 * Evaluates the collision frequency at radial grid point ir and momentum p,
 * neglecting any contribution from the nonlinear collision operator
 */
/*
real_t SlowingDownFrequency::evaluateAtP(len_t ir, real_t p){
    real_t *ncold = unknowns->GetUnknownData(id_ncold);
    real_t ntarget = ncold[ir];
    if (isNonScreened)
        ntarget += nbound[ir];

    real_t collQty;
    real_t preFact = evaluatePreFactorAtP(p); 
    collQty = lnLambdaEE->evaluateAtP(ir,p) * evaluateElectronTermAtP(ir,p) * ntarget;
    
    if(isPartiallyScreened){
        len_t ind;
        for(len_t iz = 0; iz<nZ; iz++){
            for(len_t Z0=0; Z0<Zs[iz]; Z0++){
                ind = ionIndex[iz][Z0];
                collQty +=  evaluateScreenedTermAtP(iz,Z0,p) * ionDensities[ir][ind];
            }
        }
    }
    collQty *= preFact;


    return collQty;
}
*/

real_t SlowingDownFrequency::evaluateScreenedTermAtP(len_t iz, len_t Z0, real_t p){
    len_t Z = ionHandler->GetZ(iz); 
    len_t ind = ionHandler->GetIndex(iz,Z0);
    if (atomicParameter[ind]==0)
        return 0;
    real_t gamma = sqrt(1+p*p);
    real_t beta = p/gamma;
    real_t h = p*sqrt(gamma-1)/atomicParameter[ind];
    real_t nBound = Z - Z0;
    return nBound*(log(1+pow(h,kInterpolate))/kInterpolate - beta*beta) ;

}

/**
 * Helper function to calculate a partial contribution to evaluateAtP
 */
real_t SlowingDownFrequency::evaluateElectronTermAtP(len_t ir, real_t p){
    if (collQtySettings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL){
        if(p==0)
            return ReallyLargeNumber;
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
 * Currently only for ions we have actual data for -- should investigate approximate models for other species.
 * TODO: Implement the approximate analytical 1-parameter formula (8) when data is missing, from 
 *       Sauer, Sabin, Oddershede J Chem Phys 148, 174307 (2018)
 */
real_t SlowingDownFrequency::GetAtomicParameter(len_t iz, len_t Z0){
    len_t Z = ionHandler->GetZ(iz);
    if(Z0==(Z-1)) // For hydrogenic ions the mean excitation energy is 14.9916*Z^2 eV 
        return Z*Z*14.9916 / Constants::mc2inEV;

    // Fetch value from table if it exists:
    for (len_t n=0; n<meanExcI_len; n++)
        if( (Z==meanExcI_Zs[n]) && (Z0==meanExcI_Z0s[n]) )
            return meanExcI_data[n];

    throw FVM::FVMException("Mean excitation energy for ion species: '%s' in charge state Z0 = " LEN_T_PRINTF_FMT " is missing.", ionHandler->GetName(iz).c_str(), Z0); 
}
