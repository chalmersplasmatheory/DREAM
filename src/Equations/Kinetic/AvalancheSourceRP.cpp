/**
 * Implementation of the Rosenbluth-Putvinski avalanche
 * source term, which takes the quadratic form
 *     T = S(r,p) * n_tot(r,t) * n_re(r,t),
 * where S is only a function of phase-space coordinates.
 */

#include "DREAM/Equations/Kinetic/AvalancheSourceRP.hpp"
#include "DREAM/Constants.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
AvalancheSourceRP::AvalancheSourceRP(
    FVM::Grid *kineticGrid, FVM::UnknownQuantityHandler *u,
    real_t pCutoff, real_t scaleFactor, RPSourceMode sm
) : FluidSourceTerm(kineticGrid, u), scaleFactor(scaleFactor), sourceMode(sm)
{
    id_ntot = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    id_Efield = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    this->pCutoff = pCutoff;
    // non-trivial temperature jacobian for Maxwellian-shaped particle source
    AddUnknownForJacobian(id_ntot);
    real_t e = Constants::ec;
    real_t epsmc = 4*M_PI*Constants::eps0 * Constants::me * Constants::c;
    this->preFactor = (e*e*e*e)/(epsmc*epsmc*Constants::c);
}

/**
 * Evaluates the constant (only grid dependent) source-shape function S(r,p)
 */
real_t AvalancheSourceRP::EvaluateRPSource(len_t ir, len_t i, len_t j){
    if(sourceMode == RP_SOURCE_MODE_FLUID)
        return scaleFactor*EvaluateNormalizedTotalKnockOnNumber(grid->GetRadialGrid()->GetFSA_B(ir), pCutoff);

    real_t pm = grid->GetMomentumGrid(ir)->GetP1_f(i);
    real_t pp = grid->GetMomentumGrid(ir)->GetP1_f(i+1);
    real_t dp = pp-pm;
    
    // if pCutoff lies above this cell, return 0.
    // if pCutoff lies inside this cell, set pm to pCutoff.
    if(pp<=pCutoff)
        return 0;
    else if(pm<pCutoff)
        pm = pCutoff;
     
    real_t gp = sqrt(1+pp*pp);
    real_t gm = sqrt(1+pm*pm);
    real_t pPart = ( 1/(gm-1) - 1/(gp-1) ) / dp;
    
    const real_t E = unknowns->GetUnknownData(id_Efield)[ir];
    int_t RESign = (E>=0) ? 1: -1;
    const real_t deltaHat = grid->GetAvalancheDeltaHat(ir,i,j, RESign);
    return scaleFactor*preFactor * pPart * deltaHat;
}

/**
 * Returns the source at grid point (ir,i,j).
 */
real_t AvalancheSourceRP::GetSourceFunction(len_t ir, len_t i, len_t j){
    real_t S = EvaluateRPSource(ir,i,j);
    const real_t ntot = unknowns->GetUnknownData(id_ntot)[ir];
    return S * ntot;
}

/**
 * Returns the source function at (ir,i,j) differentiated with respect to the unknown x_derivId at (ir,i,j)
 */
real_t AvalancheSourceRP::GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId){
    if(derivId==id_ntot)
        return EvaluateRPSource(ir,i,j);
    else
        return 0;
}

/**
 * Returns the flux-surface averaged avalanche source integrated over 
 * all xi and momenta pLower < p < pUpper, normalized to n_re*n_tot. 
 *  ir: radial grid index
 *  FSA_B: the flux surface average <B/Bmin> at ir
 */
real_t AvalancheSourceRP::EvaluateNormalizedTotalKnockOnNumber(real_t FSA_B, real_t pLower, real_t pUpper){
    if(pLower==0)
        return std::numeric_limits<real_t>::infinity();
    real_t e = Constants::ec;
    real_t epsmc = 4*M_PI*Constants::eps0 * Constants::me * Constants::c;
    real_t preFactor = (e*e*e*e)/(epsmc*epsmc*Constants::c);
    
    // IOverG = 1/(gamma-1)
    real_t pLo2 = pLower*pLower;
    real_t gLower = sqrt(1+pLo2);
    real_t IOverGLo = (gLower+1)/pLo2;

    real_t IOverGUp;
    if(pUpper != std::numeric_limits<real_t>::infinity()){
        real_t pUp2 = pUpper*pUpper;
        real_t gUpper = sqrt(1+pUp2);
        IOverGUp = (gUpper+1)/pUp2;
    } else 
        IOverGUp = 0;

    return 2*M_PI*preFactor*FSA_B*(IOverGLo - IOverGUp);
}