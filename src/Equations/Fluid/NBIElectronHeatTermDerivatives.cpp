/**
 * Implementation of the Jacobian derivatives for the NBI electron heating term.
 * This file handles the calculation of partial derivatives needed for the
 * implicit time evolution scheme.
 */

#include "DREAM/Equations/Fluid/NBIElectronHeatTerm.hpp"

using namespace DREAM;
/**
 * Computes the partial derivatives of the NBI heating term with respect to:
 * - Ion density (ni)
 * - Ion temperature (Ti)
 * - Electron temperature (Te)
 * - Electron density (ne)
 */
real_t NBIElectronHeatTerm::Compute_dP_derivative(len_t ir, len_t derivId, FVM::UnknownQuantityHandler *unknowns)
{
    const real_t *ncold = this->unknowns->GetUnknownData(id_ncold);
    const real_t *Tcold = this->unknowns->GetUnknownData(id_Tcold);
    const real_t *ni = this->unknowns->GetUnknownData(id_ion_density);
    const real_t *Ti = this->unknowns->GetUnknownData(id_ion_temperature);

    const real_t v_beam = std::sqrt(2.0 * this->Ti_beam / this->m_i_beam);
    ADASRateInterpolator *scd = this->adas->GetSCD(Zion);
    ADASRateInterpolator *ccd = this->adas->GetCCD(Zion);

    real_t total_deriv = 0.0;

    for (len_t i_theta = 0; i_theta < this->n_beam_theta; ++i_theta)
    {
        real_t theta_B = i_theta * this->d_beam_theta;

        for (len_t i_r = 0; i_r < this->n_beam_radius; ++i_r)
        {
            real_t r_B = i_r * this->d_beam_radius;
            real_t I_s = 0.0;
            real_t dS_dx = 0.0;

            for (len_t i_s = 0; i_s < this->n_beam_s; ++i_s)
            {
                real_t s_B = this->s_start + i_s * this->d_beam_s;
                if (r_B > this->r_beam || s_B > this->s_max)
                    continue;

                int_t ir_here = CalculatePencilBeamFindFlux(s_B, r_B, theta_B);
                if (ir_here != ir || ir_here < 0 || ir_here >= this->nr)
                    continue;

                real_t ne = ncold[ir], Te = Tcold[ir];
                real_t ni_ = ni[ir], Ti_ = Ti[ir];
                real_t I_ion = scd->Eval(this->Z0, ne, Te);
                real_t I_cx = ccd->Eval(this->Z0, ni_, Ti_);
                real_t I = I_ion + I_cx;

                real_t dI = 0.0;
                if (derivId == this->id_ncold)
                    dI = scd->Eval_deriv_n(this->Z0, ne, Te);
                else if (derivId == this->id_Tcold)
                    dI = scd->Eval_deriv_T(this->Z0, ne, Te);
                else if (derivId == this->id_ion_density)
                    dI = ccd->Eval_deriv_n(this->Z0, ni_, Ti_);
                else if (derivId == this->id_ion_temperature)
                    dI = ccd->Eval_deriv_T(this->Z0, ni_, Ti_);

                real_t lambda = v_beam / (ne * I);

                real_t dLambda = 0.0;
                if (derivId == id_ncold)
                    dLambda = -v_beam / (ne * ne * I) - v_beam * dI / (ne * I * I);
                else
                    dLambda = -v_beam * dI / (ne * I * I);

                I_s += this->d_beam_s / lambda;
                dS_dx += (1 / (lambda * lambda)) * dLambda * this->d_beam_s;

                real_t S = std::exp(-I_s);
                real_t dS = -S * dS_dx;

                real_t j_B_r = this->j_B_profile->Eval(r_B)[0];
                real_t dV = r_B * this->d_beam_radius * this->d_beam_theta * this->d_beam_s;
                real_t H_prefac = j_B_r / this->I_B * dV;

                real_t result = H_prefac * ((-1.0 / (lambda * lambda)) * dLambda * S + (1.0 / lambda) * dS);
                total_deriv += result;
            }
        }
    }

    return total_deriv * (this->beamPower / this->plasmaVolume) * this->depositedFraction;
}
