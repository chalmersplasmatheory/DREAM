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
#include <cmath>

using namespace DREAM;

/**
 * Mean excitation energy atomic data for ions from 
 * Sauer, Sabin, Oddershede J Chem Phys 148, 174307 (2018)
 */

// List of mean excitation energies in units of eV
const real_t SlowingDownFrequency::MEAN_EXCITATION_ENERGY_DATA[MAX_Z][MAX_Z] = {
/* H  */ { 14.99, NAN,   NAN,   NAN,   NAN,   NAN,   NAN,   NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* He */ { 42.68, 59.88, NAN,   NAN,   NAN,   NAN,   NAN,   NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* Li */ { 33.1,  108.3, 134.5, NAN,   NAN,   NAN,   NAN,   NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* Be */ { 42.2,  76.9,  205.0, 240.2, NAN,   NAN,   NAN,   NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* B  */ { 52.6,  82.3,  136.9, 330.4, 374.6, NAN,   NAN,   NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* C  */ { 65.9,  92.6,  134.8, 214.2, 486.2, 539.5, NAN,   NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* N  */ { 81.6,  107.4, 142.4, 200.2, 308.7, 672.0, 734.3, NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* O  */ { 97.9,  125.2, 157.2, 202.2, 278.6, 420.7, 887.8, 959.0,  NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* F  */ { 116.5, 144.0, 176.4, 215.6, 272.3, 370.2, 550.0, 1133.5, 1213.7, NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* Ne */ { 137.2, 165.2, 196.9, 235.2, 282.8, 352.6, 475.0, 696.8,  1409.2, 1498.4, NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* Na */ { 125.7, 189.2, 220.4, 256.8, 301.9, 358.7, 443.5, 593.3,  861.2,  1715.6, 1813.9, NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* Mg */ { 128.0, 173.7, 246.8, 282.5, 324.3, 376.7, 443.8, 544.8,  724.8,  1043.2, 2051.5, 2158.8, NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* Al */ { 132.2, 172.7, 225.8, 310.8, 351.0, 398.8, 459.2, 537.4,  656.4,  869.6,  1242.7, 2417.2, 2533.5, NAN,    NAN,    NAN,    NAN,    NAN},
/* Si */ { 140.8, 177.2, 221.2, 283.1, 381.4, 426.5, 480.6, 549.7,  640.1,  778.6,  1027.9, 1459.8, 2813.0, 2938.3, NAN,    NAN,    NAN,    NAN},
/* P  */ { 151.6, 185.3, 225.2, 274.3, 345.9, 458.5, 508.8, 569.7,  648.2,  751.7,  911.2,  1199.7, 1694.6, 3238.8, 3373.1, NAN,    NAN,    NAN},
/* S  */ { 162.4, 195.7, 232.8, 277.3, 332.4, 414.4, 542.1, 598.0,  666.2,  754.6,  872.2,  1054.5, 1384.9, 1947.0, 3694.5, 3837.8, NAN,    NAN},
/* Cl */ { 174.9, 206.8, 242.9, 284.1, 333.8, 395.5, 488.6, 632.1,  694.0,  769.9,  869.1,  1001.8, 1208.2, 1583.7, 2217.2, 4180.2, 4332.5, NAN},
/* Ar */ { 188.7, 219.5, 254.0, 293.7, 339.4, 394.9, 463.9, 568.6,  728.8,  797.0,  881.1,  991.6,  1140.3, 1372.6, 1796.0, 2505.0, 4695.9, 4857.2 }};

const real_t SlowingDownFrequency::MEAN_EXCITATION_ENERGY_FUNCTION_D[MAX_NE] = {0, 0.00, 0.24, 0.34, 0.41, 0.45, 0.48, 0.50, 0.51, 0.52, 0.55, 0.57, 0.58, 0.59};
const real_t SlowingDownFrequency::MEAN_EXCITATION_ENERGY_FUNCTION_S_0[MAX_NE] = {0, 0.30, 1.51, 2.32, 3.13, 3.90, 4.67, 5.44, 6.21, 6.97, 8.10, 9.08, 10.03, 10.94};
// according to Berger et al., J of the ICRU os19 22 (1984), all neutral ions with Z >= 19 have I~10*Z eV (8.8 to 11.1 eV) 
const real_t SlowingDownFrequency::HIGH_Z_EXCITATION_ENERGY_PER_Z = 10.0; 
const real_t SlowingDownFrequency::HYDROGEN_MEAN_EXCITATION_ENERGY = 14.99; // Mean excitation energy for neutral H

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
    if (Z==Z0)
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
 * When data is missing, it uses the approximate analytical 2-parameter formula (8) from 
 * Sauer, Sabin, Oddershede J Chem Phys 148, 174307 (2018)
 */
real_t SlowingDownFrequency::GetAtomicParameter(len_t iz, len_t Z0){
    len_t Z = ionHandler->GetZ(iz);

    real_t I;
    real_t D_N;
    real_t S_N0;

    if (Z == Z0){
        return NAN;
    }
    if (Z <= MAX_Z){ /* use tabulated data */
        I = MEAN_EXCITATION_ENERGY_DATA[Z-1][Z0];
    }else{ /* use the formula instead */
        len_t Ne = Z-Z0;
        if (Ne <= MAX_NE){
            D_N = MEAN_EXCITATION_ENERGY_FUNCTION_D[Ne-1]; 
            S_N0 = MEAN_EXCITATION_ENERGY_FUNCTION_S_0[Ne-1];
        }else{
            D_N = MEAN_EXCITATION_ENERGY_FUNCTION_D[MAX_NE-1]; 
            S_N0 = Ne - sqrt(Ne*HIGH_Z_EXCITATION_ENERGY_PER_Z / HYDROGEN_MEAN_EXCITATION_ENERGY); // S_N0: for a neutral atom with Z=N
        }
        real_t A_N = (1-D_N) * (1-D_N);
        real_t B_N = 2*(1-D_N) * (Ne*D_N - S_N0);
        real_t C_N = (Ne*D_N - S_N0) * (Ne*D_N - S_N0);

        I = HYDROGEN_MEAN_EXCITATION_ENERGY * (A_N*Z*Z + B_N*Z + C_N);
    }
    return I / Constants::mc2inEV;
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
    if(isPartiallyScreened)
        for(len_t iz = 0; iz<nZ; iz++)
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                len_t ind = ionIndex[iz][Z0];
                p3nuS0 += evaluateScreenedTermAtP(iz,Z0,0,collQtySettings->collfreq_mode) * ionDensities[ir][ind];
            }
    p3nuS0 *= preFactor;
    return p3nuS0;
}


/**
 * Evaluates partial derivatives of lim_{p\to 0} p^3nu_s.
 */
real_t* SlowingDownFrequency::GetPartialP3NuSAtZero(len_t derivId){
    real_t preFactor = constPreFactor;
    len_t nMultiples = 1;
    if(derivId == id_ni)
        nMultiples = nzs;
    real_t *dP3nuS = new real_t[nr*nMultiples];
    for(len_t i = 0; i<nr*nMultiples; i++)
        dP3nuS[i] = 0;

    // Set partial n_cold 
    if(derivId == id_ncold)
        for(len_t ir=0; ir<nr; ir++){
            real_t lnLee0 = lnLambdaEE->evaluateAtP(ir,0);
            dP3nuS[ir] = preFactor * lnLee0 * evaluateElectronTermAtP(ir,0,collQtySettings->collfreq_mode);
        }
    else if(derivId == id_ni)    
        for(len_t ir=0; ir<nr; ir++){
            real_t electronTerm = preFactor*evaluateElectronTermAtP(ir,0,collQtySettings->collfreq_mode);
            real_t ntarget = unknowns->GetUnknownData(id_ncold)[ir];
            if (isNonScreened)
                ntarget += ionHandler->evaluateBoundElectronDensityFromQuasiNeutrality(ir);
            for(len_t iz=0; iz<nZ; iz++)
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                    len_t indZ = ionIndex[iz][Z0];
                    dP3nuS[indZ*nr + ir] += ntarget*electronTerm*lnLambdaEE->evaluatePartialAtP(ir,0,id_ni,indZ);
                }
            if(isNonScreened)
                for(len_t iz=0; iz<nZ; iz++)
                    for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                        len_t indZ = ionIndex[iz][Z0];
                        dP3nuS[indZ*nr + ir] += (Zs[iz] - Z0) * electronTerm * lnLambdaEE->evaluateAtP(ir,0);
                    }
            else if(isPartiallyScreened)
                for(len_t iz=0; iz<nZ; iz++)
                    for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                        len_t indZ = ionIndex[iz][Z0];
                        dP3nuS[indZ*nr + ir] += preFactor * evaluateScreenedTermAtP(iz,Z0,0,collQtySettings->collfreq_mode);
                    }
        }
    else if(derivId == id_Tcold) 
        for(len_t ir=0; ir<nr; ir++){
            real_t ntarget = unknowns->GetUnknownData(id_ncold)[ir];
            if (isNonScreened)
                ntarget += ionHandler->evaluateBoundElectronDensityFromQuasiNeutrality(ir);
           
            real_t lnLee0 = lnLambdaEE->evaluateAtP(ir,0);
            real_t dLnLee0 = lnLambdaEE->evaluatePartialAtP(ir,0,id_Tcold,0);
            dP3nuS[ir] = preFactor * ntarget * (lnLee0 * evaluateDDTElectronTermAtP(ir,0,collQtySettings->collfreq_mode)
                                                + dLnLee0 * evaluateElectronTermAtP(ir,0,collQtySettings->collfreq_mode));
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
