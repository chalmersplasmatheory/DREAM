/**
 * Implementations of template specializations in the 
 * 'SvenssonTransport' equation term.
 */

#include "DREAM/Equations/SvenssonTransport.hpp"


// Template specializations
/**
 * Get advection coefficient array at the specified radius.
 */
template<>
const real_t *DREAM::SvenssonTransport<DREAM::FVM::AdvectionTerm>::GetCoefficient(
    const len_t ir
) { return this->fr[ir]; }

/**
 * Get diffusion coefficient array at the specified radius.
 */
template<>
const real_t *DREAM::SvenssonTransport<DREAM::FVM::DiffusionTerm>::GetCoefficient(
    const len_t ir
) { return this->drr[ir]; }

template<>
void DREAM::SvenssonTransport<DREAM::FVM::AdvectionTerm>::_setcoeff(
    const len_t ir, const real_t v
) { this->fr[ir][0] += v; }

template<>
void DREAM::SvenssonTransport<DREAM::FVM::DiffusionTerm>::_setcoeff(
    const len_t ir, const real_t v)
{ this->drr[ir][0] += v; }


