/**
 * Implementations of template specializations in the 
 * 'SvenssonTransport' equation term.
 */

#include "DREAM/Equations/Fluid/SvenssonTransport.hpp"

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




const real_t * SvenssonTransportDiffusionTerm::EvaluateIntegrand(len_t ir){
    // Inverse of p-bar on the Flux grid
    real_t pBarInv_f = GetPBarInv_f(ir);

    
    // Calculating the integrand
    for( len_t i=0; i < this->grid->GetNp1(0); i++){
        this->integrand[i] = -this->coeffD[ir][i] * pBarInv_f
            * exp(-(this->p[i] - this->pStar) * pBarInv_f);
    }
    return this->integrand; // An array (in p) with the value of the integrands
}


const real_t * SvenssonTransportAdvectionTermA::EvaluateIntegrand(len_t ir){
    // Inverse of p-bar on the Flux grid
    real_t pBarInv_f = GetPBarInv_f(ir);
    
    // Calculating the integrand
    for( len_t i=0; i < this->grid->GetNp1(0); i++){
        this->integrand[i] = this->coeffA[ir][i] * pBarInv_f
            * exp(-(this->p[i] - this->pStar) * pBarInv_f);
    }
    return this->integrand; // An array (in p) with the value of the integrands
}


const real_t * SvenssonTransportAdvectionTermD::EvaluateIntegrand(len_t ir){

    real_t pBarInv_f, dr_pBarInv_f;//, tmp_pBarInv_f; // Inverse of p-bar on the Flux grid

    pBarInv_f = GetPBarInv_f(ir, &dr_pBarInv_f);

    // Calculating the integrand
    for( len_t i=0; i < this->grid->GetNp1(0); i++){
        this->integrand[i] = -this->coeffD[ir][i] * pBarInv_f
            * pBarInv_f * dr_pBarInv_f * (1 - this->p[i]*pBarInv_f)
            * exp(-(this->p[i] - this->pStar) * pBarInv_f);
    }
    return this->integrand; // An array (in p) with the value of the integrands
}
