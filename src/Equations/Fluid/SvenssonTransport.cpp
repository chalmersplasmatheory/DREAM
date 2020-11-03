/**
 * Implementations of template specializations in the 
 * 'SvenssonTransport' equation term.
 */

#include "DREAM/Equations/SvenssonTransport.hpp"

using namespace DREAM;

// Template specializations
/**
 * Get advection coefficient array at the specified radius.
 */
template<>
const real_t *SvenssonTransport<FVM::AdvectionTerm>::GetCoefficient(
    const len_t ir)
{ return this->fr[ir]; }

/**
 * Get diffusion coefficient array at the specified radius.
 */
template<>
const real_t *SvenssonTransport<FVM::DiffusionTerm>::GetCoefficient(
    const len_t ir )
{ return this->drr[ir]; }

template<>
void SvenssonTransport<FVM::AdvectionTerm>::_setcoeff(
    const len_t ir, const real_t v )
{ this->fr[ir][0] += v; }

template<>
void SvenssonTransport<FVM::DiffusionTerm>::_setcoeff(
    const len_t ir, const real_t v )
{ this->drr[ir][0] += v; }




const reat_t *SvenssonTransportDiffusive::EvaluateIntegrand(len_t ir){
    // Essential values taken on the raidal grid
    real_t *E = this->unknowns->GetUnknownData(this->EID);
    real_t *EcEff = this->REFluid->effectiveCriticalField;
    real_t *tauRel = this->REFluid->tauEERel;
    real_t *gamma_r = this->REFluid->avalancheGrowthRate;

    real_t *integrand = new real_t[this->np];

    
    real_t pBarInv_f; // Inverse of p-bar on the Flux grid
    // YYY Should each individual value be interpolated or pBar itsellf?

    // Interpolating (extrapolating) the inverse of p bar onto the
    // flux grid.
    if (ir == 0) {
        // Zero flux at r = 0
        pBarInv_f = tauRel[0] * gamma_r[0] / (E[0]-EcEff[0]);
    }
    else if (ir == this->nr) {
        // Linearly extrapolating the value at the end point from the
        // two previous points.
        //
        // N.B.! The extrapolation assume that the grid cell size is
        // uniform, and that the extrapolated value lies half a grid
        // cell away from the last point.
        pBarInv_f  = 1.5 * (tauRel[ir-1] * gamma_r[ir-1] / (E[ir-1]-EcEff[ir-1]));
        pBarInv_f -= 0.5 * (tauRel[ir-2] * gamma_r[ir-2] / (E[ir-2]-EcEff[ir-2]));
    }
    else {
        // In the middle, we simple linearly interpolate
        pBarInv_f  = tauRel[ir] * gamma_r[ir] / (E[ir]-EcEff[ir]);
        pBarInv_f += tauRel[ir-1] * gamma_r[ir-1] / (E[ir-1]-EcEff[ir-1]);
        pBarInv_f *= 0.5;
    }
    
    // Calculating the integrand
    for( len_t i=0; i < this->np; i++){
        integrand[i] = -this->coeffD[ir][i] * pBarInv
            * exp(-(this->p[i] - this->pStar) * pBarInv);
    }
    return integrand; // An array (in p) with the value of the integrands
}


const reat_t *SvenssonTransportAdvectiveA::EvaluateIntegrand(len_t ir){
    // Essential values taken on the raidal grid
    real_t *E = this->unknowns->GetUnknownData(this->EID);
    real_t *EcEff = this->REFluid->effectiveCriticalField;
    real_t *tauRel = this->REFluid->tauEERel;
    real_t *gamma_r = this->REFluid->avalancheGrowthRate;

    real_t *integrand = new real_t[this->np];

    
    real_t pBarInv_f; // Inverse of p-bar on the Flux grid
    // YYY Should each individual value be interpolated or pBar itsellf?

    // Interpolating (extrapolating) the inverse of p bar onto the
    // flux grid.
    if (ir == 0) {
        // Zero flux at r = 0
        pBarInv_f = tauRel[0] * gamma_r[0] / (E[0]-EcEff[0]);
    }
    else if (ir == this->nr) {
        // Linearly extrapolating the value at the end point from the
        // two previous points.
        //
        // N.B.! The extrapolation assume that the grid cell size is
        // uniform, and that the extrapolated value lies half a grid
        // cell away from the last point.
        pBarInv_f  = 1.5 * (tauRel[ir-1] * gamma_r[ir-1] / (E[ir-1]-EcEff[ir-1]));
        pBarInv_f -= 0.5 * (tauRel[ir-2] * gamma_r[ir-2] / (E[ir-2]-EcEff[ir-2]));
    }
    else {
        // In the middle, we simple linearly interpolate
        pBarInv_f  = tauRel[ir] * gamma_r[ir] / (E[ir]-EcEff[ir]);
        pBarInv_f += tauRel[ir-1] * gamma_r[ir-1] / (E[ir-1]-EcEff[ir-1]);
        pBarInv_f *= 0.5;
    }
    
    // Calculating the integrand
    for( len_t i=0; i < this->np; i++){
        integrand[i] = this->coeffA[ir][i] * pBarInv
            * exp(-(this->p[i] - this->pStar) * pBarInv);
    }
    return integrand; // An array (in p) with the value of the integrands
}



const reat_t *SvenssonTransportAdvectiveD::EvaluateIntegrand(len_t ir){
    real_t *E = this->unknowns->GetUnknownData(this->EID);
    real_t *EcEff = this->REFluid->effectiveCriticalField;
    real_t *tauRel = this->REFluid->tauEERel;
    real_t *gamma_r = this->REFluid->avalancheGrowthRate;

    real_t *integrand = new real_t[this->np];

    real_t dr = this->grid->GetDr(ir); 
    
    real_t pBarInv_f, dr_pBarInv_f, tmp_pBarInv_f; // Inverse of p-bar on the Flux grid
    // YYY Should each individual value be interpolated or pBar itsellf?
    
    // Interpolating (extrapolating) the inverse of p bar onto the
    // flux grid.
    if (ir == 0) {
        // Zero flux at r = 0, so threrefore choose 
        pBarInv_f = tauRel[0] * gamma_r[0] / (E[0]-EcEff[0]);
        dr_pBarInv_f = 0.0;
    }
    else if (ir == this->nr) {
        // Linearly extrapolating the value at the end point from the
        // two previous points.
        //
        // N.B.! The extrapolation assume that the grid cell size is
        // uniform, and that the extrapolated value lies half a grid
        // cell away from the last point.

        // pBarInv_f  = 1.5 * (tauRel[ir-1] * gamma_r[ir-1] / (E[ir-1]-EcEff[ir-1]));
        // pBarInv_f -= 0.5 * (tauRel[ir-2] * gamma_r[ir-2] / (E[ir-2]-EcEff[ir-2]));
        pBarInv_f     = tauRel[ir-1] * gamma_r[ir-1] / (E[ir-1]-EcEff[ir-1]);
        tmp_pBarInv_f = tauRel[ir-2] * gamma_r[ir-2] / (E[ir-2]-EcEff[ir-2]);

        // N.B.! this order of operations is important
        // Using the same derivative as for the linear extrapolation.
        dr_pBarInv_f = (pBarInv_f - tmp_pBarInv_f) / dr;
        pBarInv_f *= 1.5;
        pBarInv_f -= 0.5 * tmp_pBarInv_f;
        
    }
    else {
        // In the middle, we simple linearly interpolate
        tmp_pBarInv_f = tauRel[ir-1] * gamma_r[ir-1] / (E[ir-1]-EcEff[ir-1]);
        pBarInv_f  = tauRel[ir] * gamma_r[ir] / (E[ir]-EcEff[ir]);

        // N.B.! this order of operations is important!
        dr_pBarInv_f = (pBarInv_f - tmp_pBarInv_f) / dr; // Derivative
        pBarInv_f += tmp_pBarInv_f;
        pBarInv_f *= 0.5;
    }
    // Calculating the integrand
    for( len_t i=0; i < this->np; i++){
        integrand[i] = -this->coeffD[ir][i] * pBarInv
            * pBarInv * dr_pBarInv_f * (1 - this->p[i]*pBarInv)
            * exp(-(this->p[i] - this->pStar) * pBarInv);
    }
    return integrand; // An array (in p) with the value of the integrands
}
