/**
 * Implementation of a class which handles the calculation of the slowing-down frequency nu_s.
 * The A^p component of the collision operator is given by -p*nu_s.
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
#include "DREAM/NotImplementedException.hpp"
#include "FVM/FVMException.hpp"

using namespace DREAM;

/**
 * Mean excitation energy atomic data for ions from 
 *   Sauer, S.P., Oddershede, J. and Sabin, J.R., 2015. 
 *   The mean excitation energy of atomic ions. 
 *   In Advances in Quantum Chemistry (Vol. 71, pp. 29-40). 
 *   Academic Press.
 */

// Number of mean excitation energies that there is data for (the length of lists below)
const len_t SlowingDownFrequency::meanExcI_len = 40;

// List of atomic charge numbers Z
const real_t SlowingDownFrequency::meanExcI_Zs[meanExcI_len] = 
{
/*H */ 1,
/*He*/ 2, 2,
/*Li*/ 3, 3, 3, 
/*C */ 6, 6, 6, 6, 6, 6, 
/*Ne*/ 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 
/*Ar*/ 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18
};

// List of charge states Z0
const real_t SlowingDownFrequency::meanExcI_Z0s[meanExcI_len] = 
{
/*H */ 0, 
/*He*/ 0, 1, 
/*Li*/ 0, 1, 2, 
/*C */ 0, 1, 2, 3, 4, 5, 
/*Ne*/ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
/*Ar*/ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17
};

// List of corresponding mean excitation energies in units of mc2
const real_t SlowingDownFrequency::meanExcI_data[meanExcI_len] = 
{
/*H */ 2.9295e-5, 
/*He*/ 8.3523e-05, 1.1718e-04, 
/*Li*/ 6.4775e-05, 2.1155e-04, 2.6243e-04, 
/*C */ 1.2896e-04, 1.8121e-04, 2.6380e-04, 4.1918e-04, 9.5147e-04, 0.0011, 
/*Ne*/ 2.6849e-04, 3.2329e-04, 3.8532e-04, 4.6027e-04, 5.5342e-04, 6.9002e-04, 9.2955e-04, 0.0014, 0.0028, 0.0029, 
/*Ar*/ 3.6888e-04, 4.2935e-04, 4.9667e-04, 5.7417e-04, 6.6360e-04, 7.7202e-04, 9.0685e-04, 0.0011, 0.0014, 0.0016, 0.0017, 0.0019, 0.0022, 0.0027, 0.0035, 0.0049, 0.0092, 0.0095
};

// Mean excitation energy of neutral hydrogen: also found in meanExcI_data above
const real_t SlowingDownFrequency::HYDROGEN_MEAN_EXCITATION_ENERGY = 2.9295e-5;


/**
 * Constructor
 */
SlowingDownFrequency::SlowingDownFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                CoulombLogarithm *lnLee,CoulombLogarithm *lnLei,
                enum OptionConstants::momentumgrid_type mgtype,  struct collqty_settings *cqset)
                : CollisionFrequency(g,u,ih,lnLee,lnLei,mgtype,cqset){
    hasIonTerm = false;
    gsl_ad_w = gsl_integration_workspace_alloc(1000);
}


/**
 * Destructor.
 */
SlowingDownFrequency::~SlowingDownFrequency(){
    gsl_integration_workspace_free(gsl_ad_w);
}


/**
 * Evaluates the matched Bethe formula according to Eq (2.31) in the Hesslow paper.
 * Modification: Moved the -beta^2 contribution inside the interpolation term in order
 * to preserve positivity of the contribution.
 */
real_t SlowingDownFrequency::evaluateScreenedTermAtP(len_t iz, len_t Z0, real_t p, OptionConstants::collqty_collfreq_mode collfreq_mode){
    len_t Z = ionHandler->GetZ(iz); 
    len_t ind = ionHandler->GetIndex(iz,Z0);
    if (atomicParameter[ind]==0)
        return 0;
    real_t p2 = p*p;
    real_t gamma = sqrt(1+p2);
    real_t beta2 = p2/(1+p2);
    real_t h = (p2/sqrt(1+gamma))/atomicParameter[ind];
    real_t NBound = Z - Z0;

    if (collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL)
        return NBound*log(1+pow(h*exp(-beta2),kInterpolate))/kInterpolate ;
    else 
        return NBound*log(exp(1)+h*exp(-beta2));

// previous expression as given in the paper:
//    return NBound*(log(1+pow(h,kInterpolate))/kInterpolate-beta2) ;
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
        return Z*Z*HYDROGEN_MEAN_EXCITATION_ENERGY;

    // Fetch value from table if it exists:
    for (len_t n=0; n<meanExcI_len; n++)
        if( (Z==meanExcI_Zs[n]) && (Z0==meanExcI_Z0s[n]) )
            return meanExcI_data[n];

    throw FVM::FVMException("Mean excitation energy for ion species: '%s' in charge state Z0 = " LEN_T_PRINTF_FMT " is missing.", ionHandler->GetName(iz).c_str(), Z0); 
}


/**
 * Helper function to calculate the partial contribution to evaluateAtP from free electrons
 */
real_t SlowingDownFrequency::evaluateElectronTermAtP(len_t ir, real_t p,OptionConstants::collqty_collfreq_mode collfreq_mode){
    if (collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
        if(p==0)
            return 0;
        real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
        real_t gamma = sqrt(1+p*p);
        real_t gammaMinusOne = p*p/(gamma+1); // = gamma-1
        real_t Theta = T_cold[ir] / Constants::mc2inEV;
        
        real_t M = gamma*gamma* evaluatePsi1(ir,p) - Theta * evaluatePsi0(ir,p);
        M +=  (Theta*gamma - 1) * p * exp( -gammaMinusOne/Theta );
        M /= gamma*gamma*evaluateExp1OverThetaK(Theta,2.0);
        return  M;
    } else 
        return 1;
    
}


/**
 * Helper function for integral term in bremsstrahlung formula 
 */
real_t bremsIntegrand(real_t x, void*){
    return log(1+x)/x;
}


/**
 * Evaluates the bremsstrahlung stopping power formula. Using the non-screened 
 * formula given as (4BN) in H W Koch and J W Motz, Rev Mod Phys 31, 920 (1959).
 */
real_t SlowingDownFrequency::evaluateBremsstrahlungTermAtP(len_t iz, len_t /*Z0*/, real_t p, OptionConstants::eqterm_bremsstrahlung_mode brems_mode, OptionConstants::collqty_collfreq_type /*collfreq_type*/){
    if(brems_mode != OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_STOPPING_POWER)
        return 0;
    else if(p==0)
        return 0;

    real_t preFactor = constPreFactor * Constants::alpha / (4*M_PI);
    len_t Z = ionHandler->GetZ(iz); 
    real_t gamma = sqrt(1+p*p);
    real_t gammaMinus1OverP = p/(gamma+1);
    preFactor *= Z*Z * gammaMinus1OverP;

    // The formula from ecritpaper Eq (18)
    // return preFactor * 4*M_PI*( 0.35+0.2*log(gamma) );

    real_t integralTerm,error;
    gsl_function GSL_Func;
    
    GSL_Func.function = &(bremsIntegrand);
    GSL_Func.params = nullptr; 
    real_t epsabs=0, epsrel=3e-3;
    gsl_integration_qag(&GSL_Func,0,2*p*(gamma+p),epsabs,epsrel,gsl_ad_w->limit,QAG_KEY,gsl_ad_w,&integralTerm,&error);

    real_t logTerm = log(gamma+p);
    real_t Term1 = (4.0/3.0) * (3*gamma*gamma+1)/(gamma*p) * logTerm;
    real_t Term2 = -(8*gamma+6*p)/(3*gamma*p*p)*logTerm*logTerm - 4/3;
    real_t Term3 = 2.0/(gamma*p) * integralTerm;

    return preFactor*(Term1+Term2+Term3);
}


/**
 * Helper function to calculate a partial contribution to evaluateAtP
 */
real_t SlowingDownFrequency::evaluateDDTElectronTermAtP(len_t ir, real_t p,OptionConstants::collqty_collfreq_mode collfreq_mode){
    if ( (collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL) && p){
        real_t T_cold = unknowns->GetUnknownData(id_Tcold)[ir];
        real_t gamma = sqrt(1+p*p);
        real_t gammaMinusOne = p*p/(gamma+1); // = gamma-1
        real_t Theta = T_cold / Constants::mc2inEV;
        real_t DDTheta = 1/Constants::mc2inEV;

        real_t Psi0 = evaluatePsi0(ir,p);
        real_t Psi1 = evaluatePsi1(ir,p);
        real_t Psi2 = evaluatePsi2(ir,p);
        real_t DDTPsi0 = DDTheta / (Theta*Theta) * (Psi1-Psi0);
        real_t DDTPsi1 = DDTheta / (Theta*Theta) * (Psi2-Psi1);

        real_t Denominator = gamma*gamma*evaluateExp1OverThetaK(Theta,2.0);
        real_t DDTDenominator = DDTheta/(Theta*Theta) * (gamma*gamma*evaluateExp1OverThetaK(Theta,1.0) - (1-2*Theta) * Denominator);

        real_t Numerator = gamma*gamma* Psi1 - Theta * Psi0;
        Numerator +=  (Theta*gamma - 1) * p * exp( -gammaMinusOne/Theta );
        
        real_t DDTNumerator = gamma*gamma* DDTPsi1 - (DDTheta * Psi0 + Theta * DDTPsi0 );
        DDTNumerator +=  (gamma + gammaMinusOne/(Theta*Theta) *(Theta*gamma - 1) ) * DDTheta * p * exp( -gammaMinusOne/Theta ) ;

        return  DDTNumerator  / Denominator - Numerator*DDTDenominator /(Denominator*Denominator);
    } else 
        return 0;
    
}


/**
 * Evaluates the purely momentum dependent prefactor for nu_s
 */
real_t SlowingDownFrequency::evaluatePreFactorAtP(real_t p,OptionConstants::collqty_collfreq_mode collfreq_mode){
    if(p==0) 
        return 0; 
    else if (collfreq_mode != OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_ULTRA_RELATIVISTIC)
        return constPreFactor * (1+p*p)/(p*p*p);
    else
        return constPreFactor / p;
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
    real_t p3nuS0 = lnLee0 * evaluateElectronTermAtP(ir,0,collQtySettings->collfreq_mode) * ntarget;

    // The partially screened term below will vanish identically for the 
    // formulas given in the Hesslow paper; we keep it for completeness in
    // case the model is changed in the future.
    if(isPartiallyScreened){ 
        len_t ind;
        for(len_t iz = 0; iz<nZ; iz++){
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                ind = ionIndex[iz][Z0];
                p3nuS0 +=  evaluateScreenedTermAtP(iz,Z0,0,collQtySettings->collfreq_mode) * ionDensities[ir][ind];
            }
        }
    }
    p3nuS0 *= preFactor;
    return p3nuS0;
}


/**
 * Evaluates partial derivatives of lim_{p\to 0} p^3nu_s.
 */
real_t* SlowingDownFrequency::GetPartialP3NuSAtZero(len_t derivId){
    real_t preFactor = constPreFactor;
    real_t *dP3nuS;
    // Set partial n_cold 
    if(derivId == id_ncold){
        dP3nuS = new real_t[nr];
        for(len_t ir=0; ir<nr; ir++){
            real_t lnLee0 = lnLambdaEE->evaluateAtP(ir,0);
            dP3nuS[ir] = preFactor * lnLee0 * evaluateElectronTermAtP(ir,0,collQtySettings->collfreq_mode);
        }
    } else if(derivId == id_ni){
        dP3nuS = new real_t[nr*nzs];
        for(len_t i = 0; i<nr*nzs; i++)
            dP3nuS[i] = 0;
        for(len_t ir=0; ir<nr; ir++){
            if(isNonScreened){
                real_t electronTerm = preFactor * lnLambdaEE->evaluateAtP(ir,0) * evaluateElectronTermAtP(ir,0,collQtySettings->collfreq_mode);
                for(len_t iz=0; iz<nZ; iz++)
                    for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                        len_t indZ = ionIndex[iz][Z0];
                        dP3nuS[indZ*nr + ir] += (Zs[iz] - Z0) * electronTerm;
                    }
            } else if(isPartiallyScreened){
                for(len_t iz=0; iz<nZ; iz++)
                    for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                        len_t indZ = ionIndex[iz][Z0];
                        dP3nuS[indZ*nr + ir] = preFactor * evaluateScreenedTermAtP(iz,Z0,0,collQtySettings->collfreq_mode);
                    }
            }
        }
    } else {
        dP3nuS = new real_t[nr];
        for(len_t ir=0; ir<nr; ir++)
            dP3nuS[ir] = 0;
    }
    return dP3nuS;
}


/**
 * Calculates a Rosenbluth potential matrix defined such that when it is muliplied
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
