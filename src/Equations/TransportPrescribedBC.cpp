/**
 * Implementations of template specializations in the 
 * 'TransportPrescribed' equation term.
 */

#include "DREAM/Equations/TransportPrescribedBC.hpp"


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
real_t DREAM::TransportPrescribedBC<DREAM::FVM::AdvectionTerm>::__GetSingleElement(
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
real_t DREAM::TransportPrescribedBC<DREAM::FVM::DiffusionTerm>::__GetSingleElement(
    const real_t Drr, const real_t S_wo_Drr, const real_t dr_f
) { return Drr*S_wo_Drr/dr_f; }

