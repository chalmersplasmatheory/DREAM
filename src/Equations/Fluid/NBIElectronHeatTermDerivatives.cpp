
/**
 * Implementation of de jacobian derivatives of the equation term to the heating of cold electrons from an NBI beam 
 */

#include "DREAM/Equations/Fluid/NBIElectronHeatTerm.hpp" 

namespace DREAM{

/**
 * Method that compute the 4 partical derivatives of NBI, and return them for every flux surface calucaltion
 */
std::tuple<real_t, real_t, real_t, real_t> NBIElectronHeatTerm::Compute_dP_derivative(len_t ir, FVM::UnknownQuantityHandler *unknowns)
 {
    const real_t *ncold = unknowns->GetUnknownData(id_ncold);
    const real_t *Tcold = unknowns->GetUnknownData(id_Tcold);
    const real_t *ni    = unknowns->GetUnknownData(id_ion_density);
    const real_t *Ti    = unknowns->GetUnknownData(id_ion_temperature);

    ADASRateInterpolator *scd = this->adas->GetSCD(this->Zion);
    ADASRateInterpolator *ccd = adas->GetCCD(this->Zion);

    len_t ntheta = 100;
    len_t nphi = 100;
    real_t dtheta = 2.0 * M_PI / ntheta;
    real_t dphi = 2.0 * M_PI / nphi;

    real_t total_derivative_ni = 0.0;
    real_t total_derivative_Ti = 0.0;
    real_t total_derivative_Te = 0.0;
    real_t total_derivative_ne = 0.0;
    real_t volume_element = 0.0;

    for (len_t itheta = 0; itheta < ntheta; itheta++) {
        real_t theta = itheta * dtheta;
        for (len_t iphi = 0; iphi < nphi; iphi++) {
            real_t phi = iphi * dphi;

            real_t R = radialGrid->GetFluxSurfaceR(ir, theta);
            real_t Z = radialGrid->GetFluxSurfaceZ(ir, theta);
            real_t x = R * cos(phi);
            real_t y = R * sin(phi);
            real_t z = Z;

            real_t r_B, theta_B, s_B;
            this-> CartesianToCylindrical(x, y, z, this->P0, this->n, this->a, r_B, theta_B, s_B);

            if (s_B < 0 || s_B > s_max || r_B > this->r_beam)
                continue;

            real_t ne_here = ncold[ir];
            real_t Te_here = Tcold[ir];
            real_t ni_here = ni[ir];
            real_t Ti_here = Ti[ir];

            real_t I_cx = ccd->Eval(this->Z0, ni_here, Ti_here);
            real_t I_iz = scd->Eval(this->Z0, ne_here, Te_here);
            real_t I = I_cx + I_iz;

            real_t dI_dni = ccd->Eval_deriv_n(this->Z0, ni_here, Ti_here);
            real_t dI_dTi = ccd->Eval_deriv_T(this->Z0, ni_here, Ti_here);
            real_t dI_dne = scd->Eval_deriv_n(this->Z0, ne_here, Te_here);
            real_t dI_dTe = scd->Eval_deriv_T(this->Z0, ne_here, Te_here);

            real_t survival = this->ComputeSurvivalProbability(ir, s_B, unknowns);


            real_t dSurvival_dni = 0.0;
            real_t dSurvival_dTi = 0.0;
            real_t dSurvival_dTe = 0.0;
            real_t dSurvival_dne = 0.0;

            real_t s = 0.0;
            while (s < s_B) {
                len_t ir_s = this-> FindFluxSurfaceIndexForS(s);
                if (ir_s >= nr) break;

                real_t ne_s = ncold[ir_s];
                real_t Te_s = Tcold[ir_s];
                real_t ni_s = ni[ir_s];
                real_t Ti_s = Ti[ir_s];

                dSurvival_dni += ds * ne_s * ccd->Eval_deriv_n(this->Z0, ni_s, Ti_s);
                dSurvival_dTi += ds * ne_s * ccd->Eval_deriv_T(this->Z0, ni_s, Ti_s);
                dSurvival_dTe += ds * ne_s * scd->Eval_deriv_T(this->Z0, ne_s, Te_s);
                dSurvival_dne += ds * ne_s * scd->Eval_deriv_n(this->Z0, ne_s, Te_s);
                

                s += ds;
            }

            dSurvival_dni *= -survival;
            dSurvival_dTi *= -survival;
            dSurvival_dTe *= -survival;
            dSurvival_dne *= -survival;

            real_t J = ComputeConfigurationSpaceJacobian(ir, theta);

            // Each derivative term
            real_t dH_dni = (ne_here * dI_dni * survival) + (ne_here * I * dSurvival_dni);
            real_t dH_dTi = (ne_here * dI_dTi * survival) + (ne_here * I * dSurvival_dTi);
            real_t dH_dTe = (ne_here * dI_dTe * survival) + (ne_here * I * dSurvival_dTe);
            real_t dH_dne = (ne_here * dI_dne * survival) + (ne_here * I * dSurvival_dne) - (I * survival) - (ne_here*I *I*survival) ; 
            
            total_derivative_ni += dH_dni * J * dtheta * dphi;
            total_derivative_Ti += dH_dTi * J * dtheta * dphi;
            total_derivative_Te += dH_dTe * J * dtheta * dphi;
            total_derivative_ne += dH_dne * J * dtheta * dphi;
            volume_element += J * dtheta * dphi;
        }
    }

    
    real_t K = (this-> beamPower / this-> plasmaVolume) * (this-> j_B / this-> I_B);
    total_derivative_ni *= K / volume_element;
    total_derivative_Ti *= K / volume_element;
    total_derivative_Te *= K / volume_element;
    total_derivative_ne *= K / volume_element;
    return std::make_tuple(total_derivative_ni, total_derivative_Ti, total_derivative_Te, total_derivative_ne);
;

}



//real_t dH_dne = (ne * dI_dne * survival) - (ne * I * dSurvival_dne * survival) - (I*survival) - ( ne* I *something * I_ne_sum*survival )
}