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
    real_t pBarInv_f = GetPBarInv_f(ir);
    
    //const len_t np_grid = grid->GetNp1(0);
    const len_t offset = ir*this->np;//_grid;
    // Calculating the integrand
    for( len_t i=0; i < this->np; i++){
        this->integrand[i] = -this->coeff[i+offset] * pBarInv_f
            * exp(-(this->p[i] - this->pStar) * pBarInv_f);
    }
    //return this->integrand; // An array (in p) with the value of the integrands
}


void SvenssonTransportAdvectionTermA::EvaluateIntegrand(len_t ir){
    // Inverse of p-bar on the Flux grid
    real_t pBarInv_f = GetPBarInv_f(ir);
    

    //const len_t np_grid = grid->GetNp1(0);
    const len_t offset = ir * this->np;
    // Calculating the integrand
    for( len_t i=0; i < this->np; i++){
        this->integrand[i] = this->coeff[i+offset] * pBarInv_f
            * exp(-(this->p[i] - this->pStar) * pBarInv_f);
    }
    //return this->integrand; // An array (in p) with the value of the integrands
}


void SvenssonTransportAdvectionTermD::EvaluateIntegrand(len_t ir){
    // Inverse of p-bar on the Flux grid
    real_t pBarInv_f, dr_pBarInv_f;

    pBarInv_f = GetPBarInv_f(ir, &dr_pBarInv_f);

    // const len_t np_grid = grid->GetNp1(0);
    const len_t offset = ir * this->np;//_grid;
    // Calculating the integrand
    for( len_t i=0; i < this->np; i++){
        this->integrand[i] = -this->coeff[i+offset] * pBarInv_f
            * pBarInv_f * dr_pBarInv_f * (1 - this->p[i]*pBarInv_f)
            * exp(-(this->p[i] - this->pStar) * pBarInv_f);
    }
    //return this->integrand; // An array (in p) with the value of the integrands
}
