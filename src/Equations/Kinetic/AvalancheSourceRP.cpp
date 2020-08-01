/**
 * Implementation of the Rosenbluth-Putvinski avalanche
 * source term, which takes the quadratic form
 *     T = S(r,p) * n_tot(r,t) * n_re(r,t),
 * where S is only a function of phase-space coordinates.
 */

#include "DREAM/Equations/Kinetic/AvalancheSourceRP.hpp"
#include "DREAM/Constants.hpp"

using namespace DREAM;


AvalancheSourceRP::AvalancheSourceRP(
    FVM::Grid *kineticGrid, FVM::UnknownQuantityHandler *u
) : FluidKineticSourceTerm(kineticGrid, u)
{
    id_ntot = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
        
    // non-trivial temperature jacobian for Maxwellian-shaped particle source
    AddUnknownForJacobian(id_ntot);

}

/**
 * Evaluates the constant source-shape function S(r,p)
 */
real_t AvalancheSourceRP::EvaluateRPSource(len_t ir, len_t i, len_t j){
    // placeholder
    real_t e = Constants::ec;
    real_t epsmc = Constants::eps0 * Constants::me * Constants::c;
    real_t kappa = (e*e*e*e)/(4*M_PI*epsmc*epsmc*Constants::c);
    real_t pp = grid->GetMomentumGrid(ir)->GetP1_f(i+1);
    real_t pm = grid->GetMomentumGrid(ir)->GetP1_f(i+1);
    real_t gp = sqrt(1+pp*pp);
    real_t gm = sqrt(1+pm*pm);
    real_t pPart = ( 1/(gm-1) - 1/(gp-1) ) / (pp-pm);
    const real_t deltaHat = grid->GetAvalancheDeltaHat(ir,i,j);
    return kappa * pPart * deltaHat;
}

real_t AvalancheSourceRP::GetSourceFunction(len_t ir, len_t i, len_t j){
    real_t S = EvaluateRPSource(ir,i,j);
    const real_t ntot = unknowns->GetUnknownData(id_ntot)[ir];
    return S * ntot;
}

/**
 * Returns the source function at (ir,i,j) differentiated with respect to the unknown x_derivId at (ir,i,j)
 */
real_t AvalancheSourceRP::GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId){
    real_t dS = 0;
    if(derivId==id_ntot){
        dS = EvaluateRPSource(ir,i,j);
    }
    return dS;
}


