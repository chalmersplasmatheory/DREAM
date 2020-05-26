/**
 * Implementation of a class which handles the calculation of the slowing-down frequency nu_s.
 * The A^p component of the collision operator is given by p*nu_s.
*/

/**
 * The calculations of the electron-ion contribution are based on Eq (2.31) from
 * L Hesslow et al., Generalized collision operator for fast electrons
 * interacting with partially ionized impurities, J Plasma Phys 84 (2018).
 * The relativistic thermal ee contribution is based on the expressions given in
 * Pike & Rose, Dynamical friction in a relativistic plasma, Phys Rev E 89 (2014).
 * The non-linear contribution corresponds to the isotropic component of the
 * non-relativistic operator following Rosenbluth, Macdonald & Judd, Phys Rev (1957),
 * and is described in doc/notes/theory.pdf Appendix B.
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
                enum OptionConstants::momentumgrid_type mgtype,  struct collqty_settings *cqset)
                : CollisionFrequency(g,u,ih,lnLee,lnLei,mgtype,cqset){
    hasIonTerm = false;
}

/**
 * Destructor.
 */
SlowingDownFrequency::~SlowingDownFrequency(){
//    DeallocatePartialQuantities();
//    DeallocateCollisionQuantities();
}


/**
 * Evaluates the matched Bethe formula according to Eq (2.31) in the Hesslow paper.
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
 *  Calculates a Rosenbluth potential matrix defined such that when it is muliplied
 * by the f_hot distribution vector, yields the slowing down frequency.
 */
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
    if(Z0==Z)
        return 0;
    else if(Z0==(Z-1)) // For hydrogenic ions the mean excitation energy is 14.9916*Z^2 eV 
        return Z*Z*14.9916 / Constants::mc2inEV;

    // Fetch value from table if it exists:
    for (len_t n=0; n<meanExcI_len; n++)
        if( (Z==meanExcI_Zs[n]) && (Z0==meanExcI_Z0s[n]) )
            return meanExcI_data[n];

    throw FVM::FVMException("Mean excitation energy for ion species: '%s' in charge state Z0 = " LEN_T_PRINTF_FMT " is missing.", ionHandler->GetName(iz).c_str(), Z0); 
}


/**
 * evaluates lim_{p\to 0} p^3nu_s, for use in the evaluation of the 
 * particle flux through the p=0 boundary. It corresponds to the 
 * evaluateAtP calculation evaluated at 0, if we set preFactor to 
 * constPreFactor instead of evaluatePreFactorAtP.
 */
real_t SlowingDownFrequency::GetP3NuSAtZero(len_t ir){
    real_t *ncold = unknowns->GetUnknownData(id_ncold);
    real_t ntarget = ncold[ir];
    if (isNonScreened)
        ntarget += nbound[ir];

    real_t preFactor = constPreFactor;
    real_t lnLee0 = lnLambdaEE->evaluateAtP(ir,0);
    real_t p3nuS0 = lnLee0 * evaluateElectronTermAtP(ir,0) * ntarget;

    // The partially screened term below will vanish identically for the 
    // formulas given in the Hesslow paper; we keep it for completeness in
    // case the model is changed in the future.
    if(isPartiallyScreened){ 
        len_t ind;
        for(len_t iz = 0; iz<nZ; iz++){
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                ind = ionIndex[iz][Z0];
                p3nuS0 +=  evaluateScreenedTermAtP(iz,Z0,0) * ionDensities[ir][ind];
            }
        }
    }
    p3nuS0 *= preFactor;
    return p3nuS0;
}
