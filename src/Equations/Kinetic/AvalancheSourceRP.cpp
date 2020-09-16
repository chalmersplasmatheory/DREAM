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
        return EvaluateIntegratedRPSource(ir);

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
    
    const len_t id_jhot = unknowns->GetUnknownID(OptionConstants::UQTY_J_HOT);  
    const real_t jhot = unknowns->GetUnknownData(id_jhot)[ir];
    int_t RESign;
    if(jhot>=0)
        RESign = 1;
    else
        RESign = -1;
    const real_t deltaHat = grid->GetAvalancheDeltaHat(ir,i,j, RESign);
    return scaleFactor*preFactor * pPart * deltaHat;
}

real_t AvalancheSourceRP::EvaluateIntegratedRPSource(len_t ir){
    real_t gCut = sqrt(1+pCutoff*pCutoff);
    return scaleFactor*2*M_PI*preFactor*grid->GetRadialGrid()->GetFSA_B(ir) * 1 / (gCut-1);
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
 * Returns the flux-surface averaged avalanche source 
 * integrated over all xi and over pLower < p < pUpper.
 */
real_t AvalancheSourceRP::EvaluateTotalKnockOnNumber(len_t ir, real_t pLower, real_t pUpper){
    len_t id_nre = unknowns->GetUnknownID(OptionConstants::UQTY_N_RE);
    const real_t n_re = unknowns->GetUnknownData(id_nre)[ir];
    const real_t n_tot = unknowns->GetUnknownData(id_ntot)[ir];

    real_t gUpper = sqrt(1+pUpper*pUpper);
    real_t gLower = sqrt(1+pLower*pLower);
    return scaleFactor*2*M_PI*n_re*n_tot*preFactor*grid->GetRadialGrid()->GetFSA_B(ir)*(1/(gLower-1) - 1/(gUpper-1));
}