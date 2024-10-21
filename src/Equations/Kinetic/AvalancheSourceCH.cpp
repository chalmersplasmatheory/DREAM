/**
 * Implementation of the Chiu-Harvey avalanche
 * source term, which takes the quadratic form
 *     T = S(r,p) * n_tot(r,t) * f_re(r,t),
 * where S is only a function of phase-space coordinates.
 */

#include "DREAM/Equations/Kinetic/AvalancheSourceCH.hpp"
#include "DREAM/Constants.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
AvalancheSourceCH::AvalancheSourceCH(
    FVM::Grid *kineticGrid, FVM::UnknownQuantityHandler *u,
    real_t pCutoff, real_t scaleFactor, CHSourceMode sm,
    CHSourcePitchMode sxm
) : FluidSourceTerm(kineticGrid, u), scaleFactor(scaleFactor), sourceMode(sm),
    sourceXiMode(sxm)
{
    SetName("AvalancheSourceCH");

    id_ntot = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    id_Efield = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    this->pCutoff = pCutoff;
    // non-trivial temperature jacobian for Maxwellian-shaped particle source
    AddUnknownForJacobian(u, id_ntot);
    real_t e = Constants::ec;
    real_t epsmc = Constants::eps0 * Constants::me * Constants::c;
    this->preFactor = (e*e*e*e)/(2*M_PI*epsmc*epsmc*Constants::c);
}

/**
 * Evaluates the constant (only grid dependent) source-shape function S(r,p)
 */
real_t AvalancheSourceCH::EvaluatCHSource(len_t ir, len_t i, len_t j){
    //if(sourceMode == CH_SOURCE_MODE_FLUID)
    //    return scaleFactor*EvaluateNormalizedTotalKnockOnNumber(pCutoff);

    real_t pp = grid->GetMomentumGrid(ir)->GetP1_f(i+1);
    
    // if pCutoff lies above this cell, return 0.
    // if pCutoff lies inside this cell, set pm to pCutoff.
    if(pp<=pCutoff)
        return 0;
    else if(pm<pCutoff)
        pm = pCutoff; // TODO: Have to fix this in FluxSurfaceAverager!
         
    int_t RESign;
    if (this->sourceXiMode == CH_SOURCE_PITCH_ADAPTIVE) {
        const real_t E = unknowns->GetUnknownData(id_Efield)[ir];
        RESign = (E>=0) ? 1: -1;
    } else if (this->sourceXiMode == CH_SOURCE_PITCH_POSITIVE)
        RESign = 1;
    else
        RESign = -1;

    const real_t BA = grid->GetAvalancheCHBounceAverage(ir,i,j, RESign);
    return scaleFactor * preFactor * deltaHat;
}

/**
 * Returns the source at grid point (ir,i,j).
 */
real_t AvalancheSourceCH::GetSourceFunction(len_t ir, len_t i, len_t j){
    real_t S = EvaluateCHSource(ir,i,j);
    const real_t ntot = unknowns->GetUnknownData(id_ntot)[ir];
    return S * ntot;
}

/**
 * Returns the source function at (ir,i,j) differentiated with respect to the unknown x_derivId at (ir,i,j)
 */
real_t AvalancheSourceCH::GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId){
    if(derivId==id_ntot)
        return EvaluateCHSource(ir,i,j);
    else
        return 0;
}

/** TODO: This?!
 * Returns the flux-surface averaged avalanche source integrated over 
 * all xi and momenta pLower < p < pUpper, normalized to n_re*n_tot. 
 */
real_t AvalancheSourceCH::EvaluateNormalizedTotalKnockOnNumber(real_t pLower, real_t pUpper){
    if(pLower==0)
        return std::numeric_limits<real_t>::infinity();
    real_t e = Constants::ec;
    real_t epsmc = Constants::eps0 * Constants::me * Constants::c;
    real_t preFactor = (e*e*e*e)/(2*M_PI*epsmc*epsmc*Constants::c);
    
    // IOverG = 1/(gamma-1)
    real_t pLo2 = pLower*pLower;
    real_t gLower = sqrt(1+pLo2);
    real_t IOverGLo = (gLower+1)/pLo2;

    real_t IOverGUp;
    if(pUpper == std::numeric_limits<real_t>::infinity())
        IOverGUp = 0;
    else { 
        real_t pUp2 = pUpper*pUpper;
        real_t gUpper = sqrt(1+pUp2);
        IOverGUp = (gUpper+1)/pUp2;
    }
    return 2*M_PI*preFactor*(IOverGLo - IOverGUp);
}
