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




void SvenssonTransportDiffusionTerm::EvaluateIntegrand(len_t ir){
    // Inverse of p-bar on the Flux grid
    real_t pBarInv_f = this->GetPBarInv_f(ir);
    
    const len_t offset = ir*this->np;
    // Calculating the integrand
    for( len_t i=0; i < this->np; i++){
        real_t p_tmp = p[i];
        real_t coeff_tmp = this->coeffRP[i+offset];
        real_t exp_tmp = exp( -(p_tmp - this->pStar) * pBarInv_f );
        real_t tmp = coeff_tmp * exp_tmp * pBarInv_f ;
        this->integrand[i] = tmp;
        // this->integrand[i] = this->coeffRP[i+offset] * pBarInv_f
        //     * exp(-(this->p[i] - this->pStar) * pBarInv_f);
    }
}


void SvenssonTransportAdvectionTermA::EvaluateIntegrand(len_t ir){
    // Inverse of p-bar on the Flux grid
    real_t pBarInv_f = this->GetPBarInv_f(ir);
    
    const len_t offset = ir * this->np;
    // Calculating the integrand
    for( len_t i=0; i < this->np; i++){
        this->integrand[i] = this->coeffRP[i+offset] * pBarInv_f
            * exp( -(this->p[i] - this->pStar) * pBarInv_f );
    }
    
}


void SvenssonTransportAdvectionTermD::EvaluateIntegrand(len_t ir){
    // Inverse of p-bar on the Flux grid
    real_t pBarInv_f, dr_pBarInv_f;
    pBarInv_f = this->GetPBarInv_f(ir, &dr_pBarInv_f);
    // printf("dr_pBarInv_f = %f\n",dr_pBarInv_f); fflush(stdout); // DEBUG
    
    const len_t offset = ir * this->np;
    // Calculating the integrand
    for( len_t i=0; i < this->np; i++){
        real_t pPrime = this->p[i] - this->pStar;
        // YYY Chck this sign!!!
        this->integrand[i] = -this->coeffRP[i+offset]
            // * pBarInv_f * pBarInv_f
            * dr_pBarInv_f * (1 - pPrime*pBarInv_f)
            * exp(-pPrime * pBarInv_f);
    }
}
