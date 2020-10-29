/**
 * Implementations of template specializations in the 
 * 'TransportBC' equation term.
 */

#include "DREAM/Equations/TransportBC.hpp"


// Template specializations
/**
 * Set a single element in the function vector/matrix
 * (AdvectionTerm implementation).
 *
 * Fr:       Advection coefficient.
 * S_wo_Fr:  Flux, divided by advection coefficient.
 * f(I,J,V): Function for setting vector/matrix element.
 */
template<>
real_t DREAM::TransportBC<DREAM::FVM::AdvectionTerm>::__GetSingleElement(
    const real_t Fr, const real_t S_wo_Fr, const real_t
) {
    // We use upwind interpolation...
    if (Fr > 0)
        return Fr*S_wo_Fr;
    else
        return 0.0;
}

/**
 * Set a single element in the function vector/matrix
 * (AdvectionTerm implementation).
 *
 * Drr:      Advection coefficient.
 * S_wo_Drr: Flux, divided by advection coefficient.
 * f(I,J,V): Function for setting vector/matrix element.
 */
template<>
real_t DREAM::TransportBC<DREAM::FVM::DiffusionTerm>::__GetSingleElement(
    const real_t Drr, const real_t S_wo_Drr, const real_t dr_f
) { return Drr*S_wo_Drr/dr_f; }

/**
 * Returns either an advection or diffusion coefficient,
 * depending on the type of the object.
 */
template<>
const real_t *DREAM::TransportBC<DREAM::FVM::AdvectionTerm>::GetCoefficient(const len_t ir) {
    return this->transportOperator->GetAdvectionCoeffR(ir);
}

template<>
const real_t *DREAM::TransportBC<DREAM::FVM::DiffusionTerm>::GetCoefficient(const len_t ir) {
    return this->transportOperator->GetDiffusionCoeffRR(ir);
}

