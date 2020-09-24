/**
 * Implementations of template specializations in the 
 * 'TransportPrescribed' equation term.
 */

#include "DREAM/Equations/TransportPrescribed.hpp"


// Template specializations
/**
 * Get advection coefficient array at the specified radius.
 */
template<>
const real_t *DREAM::TransportPrescribed<DREAM::FVM::AdvectionTerm>::GetCoefficient(
    const len_t ir
) { return this->fr[ir]; }

/**
 * Get diffusion coefficient array at the specified radius.
 */
template<>
const real_t *DREAM::TransportPrescribed<DREAM::FVM::DiffusionTerm>::GetCoefficient(
    const len_t ir
) { return this->drr[ir]; }

template<>
void DREAM::TransportPrescribed<DREAM::FVM::AdvectionTerm>::_setcoeff(
    const len_t ir, const len_t j, const real_t v
) { this->fr[ir][j] = v; }

template<>
void DREAM::TransportPrescribed<DREAM::FVM::DiffusionTerm>::_setcoeff(
    const len_t ir, const len_t j, const real_t v
) { this->drr[ir][j] = v; }

