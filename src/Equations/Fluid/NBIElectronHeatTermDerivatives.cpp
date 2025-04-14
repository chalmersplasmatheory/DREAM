/**
 * Implementation of the Jacobian derivatives for the NBI electron heating term.
 * This file handles the calculation of partial derivatives needed for the 
 * implicit time evolution scheme.
 */

#include "DREAM/Equations/Fluid/NBIElectronHeatTerm.hpp" 

namespace DREAM {

/**
 * Computes the partial derivatives of the NBI heating term with respect to:
 * - Ion density (ni)
 * - Ion temperature (Ti)
 * - Electron temperature (Te)
 * - Electron density (ne)
 */
std::tuple<real_t, real_t, real_t, real_t> NBIElectronHeatTerm::Compute_dP_derivative(
    len_t ir, FVM::UnknownQuantityHandler *unknowns)
{
    // Get plasma profiles at this radius
    const real_t *ncold = unknowns->GetUnknownData(id_ncold);
    const real_t *Tcold = unknowns->GetUnknownData(id_Tcold);
    const real_t *ni    = unknowns->GetUnknownData(id_ion_density);
    const real_t *Ti    = unknowns->GetUnknownData(id_ion_temperature);
    this-> ds= ds;

    // Get rate coefficients from ADAS
    ADASRateInterpolator *scd = this->adas->GetSCD(this->Zion);  // Ionization
    ADASRateInterpolator *ccd = adas->GetCCD(this->Zion);        // Charge exchange

    // Set up angular grid for integration
    len_t ntheta = 100;
    len_t nphi = 100;
    real_t dtheta = 2.0 * M_PI / ntheta;
    real_t dphi = 2.0 * M_PI / nphi;

    // Initialize accumulators for derivatives and volume
    real_t total_derivative_ni = 0.0;
    real_t total_derivative_Ti = 0.0;
    real_t total_derivative_Te = 0.0;
    real_t total_derivative_ne = 0.0;
    real_t volume_element = 0.0;

    // Loop over poloidal and toroidal angles
    for (len_t itheta = 0; itheta < ntheta; itheta++) {
        real_t theta = itheta * dtheta;
        for (len_t iphi = 0; iphi < nphi; iphi++) {
            real_t phi = iphi * dphi;

            // Calculate position in tokamak coordinates
            real_t R = this->R0 + radialGrid->GetFluxSurfaceR(ir, theta);
            real_t Z = radialGrid->GetFluxSurfaceZ(ir, theta);
            real_t x = R * cos(phi);
            real_t y = R * sin(phi);
            real_t z = Z;

            // Transform to beam coordinates
            real_t r_B, theta_B, s_B;
            this->CartesianToCylindrical(x, y, z, this->P0, this->n, this->a, 
                                       r_B, theta_B, s_B);

            // Skip points outside beam
            if (s_B < 0 || s_B > s_max || r_B >= this->r_beam)
                continue;

            // Get local plasma parameters
            real_t ne_here = ncold[ir];
            real_t Te_here = Tcold[ir];
            real_t ni_here = ni[ir];
            real_t Ti_here = Ti[ir];

            // Calculate interaction rates and their derivatives
            real_t I_cx = ccd->Eval(this->Z0, ni_here, Ti_here);
            real_t I_iz = scd->Eval(this->Z0, ne_here, Te_here);
            real_t I = I_cx + I_iz;

            real_t dI_dni = ccd->Eval_deriv_n(this->Z0, ni_here, Ti_here);
            real_t dI_dTi = ccd->Eval_deriv_T(this->Z0, ni_here, Ti_here);
            real_t dI_dne = scd->Eval_deriv_n(this->Z0, ne_here, Te_here);
            real_t dI_dTe = scd->Eval_deriv_T(this->Z0, ne_here, Te_here);

            // Calculate beam survival probability and its derivatives
            real_t survival = this->ComputeSurvivalProbability(s_B, r_B, theta_B, unknowns);

            // Initialize survival probability derivatives
            real_t dSurvival_dni = 0.0;
            real_t dSurvival_dTi = 0.0;
            real_t dSurvival_dTe = 0.0;
            real_t dSurvival_dne = 0.0;

            // Integrate along beam path to calculate survival derivatives
            real_t s = 0.0;
            while (s < s_B) {
                len_t ir_s = this->CalculatePencilBeamFindFlux(s, r_B, theta_B);
                if (ir_s >= nr) break;

                // Get plasma parameters along beam path
                real_t ne_s = ncold[ir_s];
                real_t Te_s = Tcold[ir_s];
                real_t ni_s = ni[ir_s];
                real_t Ti_s = Ti[ir_s];

                // Accumulate derivatives of survival probability
                dSurvival_dni += ds * ne_s * ccd->Eval_deriv_n(this->Z0, ni_s, Ti_s);
                dSurvival_dTi += ds * ne_s * ccd->Eval_deriv_T(this->Z0, ni_s, Ti_s);
                dSurvival_dTe += ds * ne_s * scd->Eval_deriv_T(this->Z0, ne_s, Te_s);
                dSurvival_dne += ds * ne_s * scd->Eval_deriv_n(this->Z0, ne_s, Te_s);

                s += ds;
            }

            // Apply chain rule to survival probability derivatives
            dSurvival_dni *= -survival;
            dSurvival_dTi *= -survival;
            dSurvival_dTe *= -survival;
            dSurvival_dne *= -survival;

            // Get Jacobian for configuration space integration
            real_t J = this->radialGrid->ComputeConfigurationSpaceJacobian(ir, theta);

            // Calculate total derivatives of heating term
            real_t dH_dni = (ne_here * dI_dni * survival) + (ne_here * I * dSurvival_dni);
            real_t dH_dTi = (ne_here * dI_dTi * survival) + (ne_here * I * dSurvival_dTi);
            real_t dH_dTe = (ne_here * dI_dTe * survival) + (ne_here * I * dSurvival_dTe);
            real_t dH_dne = (ne_here * dI_dne * survival) + (ne_here * I * dSurvival_dne) 
                           - (I * survival) - (ne_here * I * I * survival);

            // Accumulate derivatives with geometric factors
            total_derivative_ni += dH_dni * J * dtheta * dphi;
            total_derivative_Ti += dH_dTi * J * dtheta * dphi;
            total_derivative_Te += dH_dTe * J * dtheta * dphi;
            total_derivative_ne += dH_dne * J * dtheta * dphi;
            volume_element += J * dtheta * dphi;
        }
    }

    // Apply normalization factors
    const real_t* j_B_values = this->j_B_profile->Eval(ir);
    real_t K = (this->beamPower / this->plasmaVolume) * (j_B_values[0] / this->I_B);
    total_derivative_ni *= K / volume_element;
    total_derivative_Ti *= K / volume_element;
    total_derivative_Te *= K / volume_element;
    total_derivative_ne *= K / volume_element;

    return std::make_tuple(total_derivative_ni, total_derivative_Ti, 
                          total_derivative_Te, total_derivative_ne);
}

} // namespace DREAM