/**
 * Implementation of equation term to the heating of cold electrons from an NBI beam
 */
#include <cmath>
#include <set>
#include <array>
#include "DREAM/Equations/Fluid/NBIElectronHeatTerm.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "DREAM/ADAS.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Equations/Fluid/IonRateEquation.hpp"
#include <tuple>
#include <algorithm>
#include "FVM/UnknownQuantityHandler.hpp"
#include <omp.h>
#include <unistd.h>
#include "FVM/Interpolator1D.hpp"
#include "DREAM/DREAMException.hpp"
#include <iostream>

using namespace DREAM;
/**
 * Constructor
 */
NBIElectronHeatTerm::NBIElectronHeatTerm(
    FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns, ADAS *adas, IonHandler *ions,
    real_t s_max, real_t r_beam, const real_t P0[3], const real_t n[3],
    real_t Ti_beam, real_t m_i_beam, real_t beamPower, FVM::Interpolator1D *j_B_profile, 
    real_t Z0, real_t Zion, real_t R0, bool TCVGaussian, FVM::Interpolator1D *Power_Profile
) : EquationTerm(grid){

    // Set grid
    this->radialGrid = grid->GetRadialGrid();
    this->nr = grid->GetNr();
    this->NBIHeatTerm = new real_t[nr];
    this->Deposition_profile = new real_t[nr];
    this->Deposition_profile_times_Vprime = new real_t[nr];
    this->H_r_dTe = new real_t[nr];
    this->H_r_dni = new real_t[nr];
    this->H_r_dTi = new real_t[nr];
    this->H_r_dne = new real_t[nr];
    this->dV_beam_prime_tot = new real_t[nr];

    // Get from other classes
    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_ion_density = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    this->id_ion_temperature = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->nZ = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this->adas = adas;
    this->ions = ions;

    // Tell DREAM these are dependencies
    this->AddUnknownForJacobian(unknowns, id_ncold);
    this->AddUnknownForJacobian(unknowns, id_Tcold);
    this->AddUnknownForJacobian(unknowns, id_ion_density);
    this->AddUnknownForJacobian(unknowns, id_ion_temperature);

    // Set integration step size
    this->unknowns = unknowns;
    this->ntheta = 50;
    this->dtheta = 2.0 * M_PI / ntheta;

    // Set beam parameters
    this->s_max = s_max;
    this->r_beam = r_beam;
    this->Ti_beam = Ti_beam;
    this->m_i_beam = m_i_beam;
    this->Z0 = Z0;
    this->Zion = Zion;
    this->R0 = R0;
    if (this->P0 == nullptr || this->n == nullptr)
            throw DREAMException("Either direction or starting point of the beam is not set correctly");
    for (int i = 0; i < 3; ++i){
        this->P0[i] = P0[i];
        this->n[i] = n[i];
    }
    

    this->beamPower = beamPower;

    // Compute the beam basis vectors e1 and e2 and the start and stop of the beam in beam coordinates, initialized to extreme values
    this->s_start = 100;
    this->s_stop = 0;
    PrecomputeBeamBasisVectors();
    PrecomputeFluxSurfaces();
    if (s_start > s_stop){
        throw DREAMException("Beam does not intersect with any flux surface.");
    }
    if (s_max < s_stop){
        throw DREAMException("Beam does not extend all the way though the plasma.");
    }

    this->j_B_profile = j_B_profile;
    this->Power_Profile = Power_Profile;
    this->n_beam_radius = 25;
    this->d_beam_radius = r_beam / n_beam_radius;
    this->n_beam_theta = 25;
    this->d_beam_theta = 2.0 * M_PI / n_beam_theta;
    this->n_beam_s = 50;
    this->d_beam_s = (s_stop - s_start) / n_beam_s;
    this->TCVGaussian = TCVGaussian;

    // Integrate VpVol to get plasma volume
    this->plasmaVolume = 0;
    const real_t *VpVol = grid->GetVpVol();
    const real_t *dr = this->radialGrid->GetDr();
    for (len_t ir = 0; ir < nr; ir++)
        this->plasmaVolume += VpVol[ir] * dr[ir] * R0;
    this->beamVolume = (s_stop - s_start) * M_PI * r_beam * r_beam;

    // Trapezoidal integration to get total beam current
    this->I_B = 0.0;
    const real_t *r_beam_grid = this->j_B_profile->GetX();
    len_t N = this->j_B_profile->GetNx();
    for (len_t i = 0; i < N - 1; ++i){
        real_t r1 = r_beam_grid[i];
        real_t r2 = r_beam_grid[i + 1];
        real_t dr = r2 - r1;

        // Evaluate current density at physical r values
        real_t j1 = this->j_B_profile->Eval(r1)[0];
        real_t j2 = this->j_B_profile->Eval(r2)[0];

        // Trapezoidal integration
        this->I_B += 2.0 * M_PI * 0.5 * (j1 * r1 + j2 * r2) * dr;
    }
}

/**
 * Calculate the beam current density j_B divided by I_B.
 */
real_t NBIElectronHeatTerm::Calculate_jB_IB(real_t r_B, real_t theta_B){
    if (TCVGaussian){
        real_t Dx = 0.1;
        real_t Dy = 0.094;
        real_t A = 4 / (M_PI * Dx * Dy);
        real_t cs = std::cos(theta_B);
        real_t sn = std::sin(theta_B);
        real_t jB_divided_IB = A * std::exp(-4*r_B*r_B * (cs*cs/(Dx*Dx) + sn*sn/(Dy*Dy)));
        return jB_divided_IB;
    } else {
        const real_t *j_B_values = j_B_profile->Eval(r_B);
        real_t j_B_r = j_B_values[0];
        return j_B_r / this->I_B;
    }
}

/**
 * Precomputes the flux surfaces in R and Z. Also find the start and end of the beam in beam coordinates.
 */
void NBIElectronHeatTerm::PrecomputeFluxSurfaces(){

    cachedFluxSurfaces.resize(nr);
    for (len_t ir = 0; ir < nr; ++ir){
        cachedFluxSurfaces[ir].resize(ntheta);
        for (real_t itheta = 0; itheta < ntheta; ++itheta){
            real_t theta = itheta * dtheta;
            real_t R = R0 + radialGrid->GetFluxSurfaceRMinusR0_theta(ir, theta);
            real_t Z = radialGrid->GetFluxSurfaceZMinusZ0_theta(ir, theta);
            cachedFluxSurfaces[ir][itheta] = {R, Z};

            real_t nphi = 100;
            real_t dphi = 2.0 * M_PI / nphi;
            if (ir == nr - 1){ // Most outermost surface
                for (real_t iphi = 0; iphi < nphi; ++iphi){
                    real_t phi = iphi * dphi;
                    real_t x = R * std::cos(phi);
                    real_t y = R * std::sin(phi);
                    real_t z = Z;
                    real_t r_B, theta_B, s_B;
                    CartesianToCylindrical(x, y, z, this->P0, this->n, r_B, theta_B, s_B);
                    if (r_B < r_beam){
                        // Save the largest value as s_stop and the lowest as s_start
                        if (s_B > s_stop)
                            s_stop = s_B;
                        if (s_B < s_start)
                            s_start = s_B;
                    }
                }
            }
        }
    }
}

/**
 * Precomputes the beam basis vectors of beam e1 and e2.
 */
void NBIElectronHeatTerm::PrecomputeBeamBasisVectors(){
    real_t norm_n = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (int i = 0; i < 3; ++i)
        n_norm[i] = n[i] / norm_n;

    std::array<real_t, 3> a;
    if (std::abs(n_norm[2]) < 0.9)
        a = {0, 0, 1};
    else
        a = {1, 0, 0};
    real_t a_dot_n = a[0] * n_norm[0] + a[1] * n_norm[1] + a[2] * n_norm[2];
    for (int i = 0; i < 3; ++i)
        e1[i] = a[i] - a_dot_n * n_norm[i];

    real_t norm_e1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
    for (int i = 0; i < 3; ++i)
        e1[i] /= norm_e1;

    e2[0] = n_norm[1] * e1[2] - n_norm[2] * e1[1];
    e2[1] = n_norm[2] * e1[0] - n_norm[0] * e1[2];
    e2[2] = n_norm[0] * e1[1] - n_norm[1] * e1[0];
}

/**
 * Destructor
 */
NBIElectronHeatTerm::~NBIElectronHeatTerm(){
    delete[] NBIHeatTerm;
    delete[] Deposition_profile;
    delete[] Deposition_profile_times_Vprime;
    delete[] H_r_dTe;
    delete[] H_r_dni;
    delete[] H_r_dTi;
    delete[] H_r_dne;
}

/**
 * Rebuild: Called at the start of every timestep.
 * Compute added heating term for each timestep and flux surface radius
 */
void NBIElectronHeatTerm::Rebuild(const real_t t, const real_t, FVM::UnknownQuantityHandler *unknowns){
    // Setting values to zero
    for (len_t ir = 0; ir < nr; ++ir){
        NBIHeatTerm[ir] = 0.0;
        Deposition_profile[ir] = 0.0;
        Deposition_profile_times_Vprime[ir] = 0.0;
        dV_beam_prime_tot[ir] = 0.0;
        H_r_dTe[ir] = 0.0;
        H_r_dni[ir] = 0.0;
        H_r_dTi[ir] = 0.0;
        H_r_dne[ir] = 0.0;
    }
    // Set the beam power, either constant or time-dependent
    real_t powerNow = 0;
    // If the time dependant value is set:
    if (this->Power_Profile != nullptr){
        const real_t *p = this->Power_Profile->Eval(t);
        powerNow = p[0];
    } else { // If the constant value is set:
        if (this->beamPower > 0)
            powerNow = this->beamPower;
        else
            throw DREAMException("No beam power set (neither constant nor profile)");
        
    }
    // If the power is not 0, then we want to compute the caluclation
    if (powerNow != 0) {
        // H(r)
        ComputeDepositionProfile(unknowns);

        // Compute the deposited fraction of the beam
        real_t cumulative = 0.0;
        real_t frac = 0.0;
        for (len_t ir = 0; ir < nr; ir++){
            cumulative += Deposition_profile_times_Vprime[ir] * d_beam_radius;
            frac = cumulative;
        }

        for (len_t ir = 0; ir < nr; ++ir){

            H_r_dTe[ir] *=  -powerNow * d_beam_radius / (radialGrid->GetVpVol()[ir]*radialGrid->GetDr()[ir]);
            H_r_dTi[ir] *= -powerNow * d_beam_radius / (radialGrid->GetVpVol()[ir]*radialGrid->GetDr()[ir]);
            H_r_dni[ir] *= -powerNow * d_beam_radius / (radialGrid->GetVpVol()[ir]*radialGrid->GetDr()[ir]);
            H_r_dne[ir] *= -powerNow * d_beam_radius / (radialGrid->GetVpVol()[ir]*radialGrid->GetDr()[ir]);
        }

        for (len_t ir = 0; ir < nr; ir++){
            NBIHeatTerm[ir] = powerNow * Deposition_profile_times_Vprime[ir] * d_beam_radius /  (radialGrid->GetVpVol()[ir]*radialGrid->GetDr()[ir]); 
            NBIHeatTerm[0] = 0;            
            
        }
    
    } else { // If power is 0
        for (len_t ir = 0; ir < nr; ir++){
            NBIHeatTerm[ir] = 0;
            H_r_dTe[ir] = 0;
            H_r_dTi[ir] = 0;
            H_r_dni[ir] = 0;
            H_r_dne[ir] = 0;
        }
    }
}

/**
 * Computes and mean free path variable for the NBI heating term
 */
void NBIElectronHeatTerm::ComputeMeanFreePath(
    len_t ir, real_t ncold, real_t Tcold, 
    real_t ni, real_t Ti, real_t &lambda_s, 
    real_t &dlambda_dI, real_t &dlambda_dne, 
    real_t &dI_dni, real_t &dI_dTi, real_t &dI_dTe
){
        // Get the SCD rate coefficient 
        ADASRateInterpolator *scd = adas->GetSCD(Zion);
        real_t dI_dne = scd->Eval_deriv_n(Z0, ncold, Tcold);
        real_t I_ion = scd->Eval(Z0, ncold, Tcold);
        
        dI_dTe = scd->Eval_deriv_T(Z0, ncold, Tcold);
        real_t v_NBI = std::sqrt(2.0 * Ti_beam / m_i_beam);
        real_t total_CX = 0.0;

        // Loop over all atomic species
        for (len_t iz = 0; iz < ions->GetNZ(); iz++){
            len_t Zmax = ions->GetZ(iz);
            
            // Use the CCD table for this nuclear charge
            ADASRateInterpolator *ccd = adas->GetCCD(Zmax);

            // Loop over charge states of this species
            for (len_t Z0 = 0; Z0 <= Zmax; Z0++){
                real_t niZ = ions->GetIonDensity(ir, iz, Z0);
                if (niZ <= 0)
                    continue;
                
                // Evaluate CCD rate coefficient
                real_t rateCoeff = ccd->Eval(Z0, niZ, Ti); //ti
                total_CX += rateCoeff* niZ; 

                dI_dni += ccd->Eval_deriv_n(Z0, niZ, Ti);
                dI_dTi= ccd->Eval_deriv_T(Z0, niZ, Ti);
            }
        }
        
        lambda_s = v_NBI / (ncold * (I_ion) + total_CX); 
        dlambda_dI = -v_NBI / (ncold * std::pow(I_ion + total_CX, 2));
        dlambda_dne = v_NBI * (-(I_ion + total_CX + dI_dne * ncold) / (ncold * ncold * (I_ion + total_CX) * (I_ion + total_CX)));
            
}

/**
 * Computes the flux surface index for a given pencil beam position (s_B, r_B, theta_B).
 * This function checks if the pencil beam intersects with any flux surface and returns the index of the first one it finds.
 */
int_t NBIElectronHeatTerm::CalculatePencilBeamFindFlux(real_t s_B, real_t r_B, real_t theta_B){
    // Convert the beam coordinates (s_B, r_B, theta_B) to Cartesian coordinates
    real_t dx = r_B * cos(theta_B);
    real_t dy = r_B * sin(theta_B);
    std::array<real_t, 3> P0_current = {
        P0[0] + dx * e1[0] + dy * e2[0],
        P0[1] + dx * e1[1] + dy * e2[1],
        P0[2] + dx * e1[2] + dy * e2[2]
    };

    int_t best_ir = -1;
    real_t smallest_diff = std::numeric_limits<real_t>::infinity();

    for (len_t ir = 0; ir < nr; ir++){ 

        for (len_t i = 0; i < ntheta; i++){ 
            len_t j = (i + 1) % ntheta;

            const FluxSurfacePoint &pt1 = cachedFluxSurfaces[ir][i];
            const FluxSurfacePoint &pt2 = cachedFluxSurfaces[ir][j];

            real_t R1 = pt1.R, Z1 = pt1.Z;
            real_t R2 = pt2.R, Z2 = pt2.Z;

            real_t dZ = Z2 - Z1;
            real_t dR = R2 - R1;

            real_t ratio = dR / dZ;

            real_t a0 = R1 * R1 - P0_current[0] * P0_current[0] - P0_current[1] * P0_current[1] + 2 * R1 * (P0_current[2] - Z1) * ratio + (P0_current[2] - Z1) * (P0_current[2] - Z1) * ratio * ratio;

            real_t a1 = 2 * R1 * n_norm[2] * ratio + 2 * n_norm[2] * (P0_current[2] - Z1) * ratio * ratio - 2 * (P0_current[0] * n_norm[0] + P0_current[1] * n_norm[1]);

            real_t a2 = n_norm[2] * n_norm[2] * ratio * ratio - n_norm[0] * n_norm[0] - n_norm[1] * n_norm[1];

            real_t discriminant = a1 * a1 / (4 * a2 * a2) - a0 / a2;
            if (discriminant <= 0)
                continue;
            

            real_t sqrt_disc = std::sqrt(discriminant);
            real_t l1 = (-a1 / (2 * a2) + sqrt_disc);
            real_t l2 = (-a1 / (2 * a2) - sqrt_disc);

            // Check both solutions
            auto check_solution = [&](real_t l) -> std::pair<bool, real_t>{
                if (l < 0)
                    return {false, 0};
                real_t t = (P0_current[2] - Z1 + l * n_norm[2]) / dZ;
                if (t < 0 || t > 1)
                    return {false, 0};
                return {true, std::abs(l - s_B)};
            };

            // Check both solutions and update if better match found
            auto [valid1, diff1] = check_solution(l1);

            if (valid1 && diff1 < smallest_diff){
                smallest_diff = diff1;
                best_ir = ir;
            }

            auto [valid2, diff2] = check_solution(l2);
            if (valid2 && diff2 < smallest_diff){
                smallest_diff = diff2;
                best_ir = ir;
            }
        }
    }
    return best_ir;
}
/**
 * Takes in the beam paramaters and turns a point on a flux surface (in x y z) into the beam coordinate system
 */
void NBIElectronHeatTerm::CartesianToCylindrical(
    real_t x, real_t y, real_t z,
    const real_t P0[3], const real_t n[3],
    real_t &r, real_t &theta, real_t &s
) {

    std::array<real_t, 3> d = {x - P0[0], y - P0[1], z - P0[2]};
    s = d[0] * n_norm[0] + d[1] * n_norm[1] + d[2] * n_norm[2];

    std::array<real_t, 3> d_trans;
    for (int i = 0; i < 3; ++i)
        d_trans[i] = d[i] - s * n_norm[i];

    real_t r_cos_theta = d_trans[0] * e1[0] + d_trans[1] * e1[1] + d_trans[2] * e1[2];
    real_t r_sin_theta = d_trans[0] * e2[0] + d_trans[1] * e2[1] + d_trans[2] * e2[2];
    r = std::sqrt(r_cos_theta * r_cos_theta + r_sin_theta * r_sin_theta);
    theta = atan2(r_sin_theta, r_cos_theta);
    if (theta < 0)
        theta += 2 * M_PI;
}

/**
 * Computes the deposition profile for a flux surface
 */
void NBIElectronHeatTerm::ComputeDepositionProfile(FVM::UnknownQuantityHandler *unknowns){

    for (len_t ir = 0; ir < nr; ++ir){
        Deposition_profile_times_Vprime[ir] = 0.0;
        Deposition_profile[ir] = 0.0;
        dV_beam_prime_tot[ir] = 0.0;
        H_r_dTe[ir] = 0.0;
        H_r_dni[ir] = 0.0;
        H_r_dTi[ir] = 0.0;
        H_r_dne[ir] = 0.0;
    }
    std::vector<real_t> lambda_cache(nr, -1.0); // -1 = not yet computed
    std::vector<real_t> dlambda_dI_cache(nr, -1.0);
    std::vector<real_t> dlambda_dne_cache(nr, -1.0);

    for (real_t i_beam_theta = 0; i_beam_theta < n_beam_theta; i_beam_theta++){
        real_t beam_theta = i_beam_theta * d_beam_theta;
        for (real_t i_beam_radius = 0; i_beam_radius < n_beam_radius; i_beam_radius++){

            real_t beam_radius = (i_beam_radius + 0.5) * d_beam_radius;
            real_t I_s = 0.0;
            real_t I_s_squared = 0.0;
            for (real_t i_beam_s = 0; i_beam_s < n_beam_s; i_beam_s++){
                real_t beam_s = s_start + i_beam_s * d_beam_s;
                len_t ir_now = CalculatePencilBeamFindFlux(beam_s, beam_radius, beam_theta);
                if (ir_now < 0 || ir_now >= nr){
                    throw DREAMException(
                        "NBIElectronHeatTerm: ir_now out of bounds: "
                        "ir_now = " LEN_T_PRINTF_FMT
                        " (s_B = %.4f, r_B = %.4f, theta_B = %.4f)",
                        ir_now, beam_s, beam_radius, beam_theta
                    );
                }
                real_t ncold = this->unknowns->GetUnknownData(id_ncold)[ir_now];
                real_t Tcold = this->unknowns->GetUnknownData(id_Tcold)[ir_now];
                real_t ni = this->unknowns->GetUnknownData(id_ion_density)[ir_now];
                real_t Ti = this->unknowns->GetUnknownData(id_ion_temperature)[ir_now];

                real_t lambda_s, dlambda_dI, dlambda_dne, dI_dni, dI_dTi, dI_dTe;
                ComputeMeanFreePath(ir_now, ncold, Tcold, ni, Ti, lambda_s, dlambda_dI, dlambda_dne, dI_dni, dI_dTi, dI_dTe);

                if (lambda_cache[ir_now] < 0){
                    lambda_cache[ir_now] = lambda_s;
                    dlambda_dI_cache[ir_now] = dlambda_dI;
                    dlambda_dne_cache[ir_now] = dlambda_dne;
                } else {
                    lambda_s = lambda_cache[ir_now];
                    dlambda_dI = dlambda_dI_cache[ir_now];
                    dlambda_dne = dlambda_dne_cache[ir_now];
                }
               
                
                I_s += 1 / lambda_s * d_beam_s;
                I_s_squared += 1 / (lambda_s * lambda_s) * d_beam_s;
                real_t survivalProb = std::exp(-I_s);

                real_t jB_divided_IB = Calculate_jB_IB(beam_radius, beam_theta);

                real_t dV_beam_prime = beam_radius * d_beam_theta * d_beam_s;
                real_t H_r_theta_phi = jB_divided_IB * (1.0 / lambda_s) * survivalProb;
                Deposition_profile[ir_now] += H_r_theta_phi;

                Deposition_profile_times_Vprime[ir_now] += H_r_theta_phi * dV_beam_prime;

                dV_beam_prime_tot[ir_now] += dV_beam_prime;

                // Derivatives
                H_r_dTe[ir_now] += jB_divided_IB * dlambda_dI * dI_dTe * ((-1.0 / (lambda_s * lambda_s)) * survivalProb + 1 / lambda_s * I_s_squared * survivalProb) * dV_beam_prime;
                H_r_dni[ir_now] += jB_divided_IB * dlambda_dI * dI_dni * ((-1.0 / (lambda_s * lambda_s)) * survivalProb + 1 / lambda_s * I_s_squared * survivalProb) * dV_beam_prime;
                H_r_dTi[ir_now] += jB_divided_IB * dlambda_dI * dI_dTi * ((-1.0 / (lambda_s * lambda_s)) * survivalProb + 1 / lambda_s * I_s_squared * survivalProb) * dV_beam_prime;
                H_r_dne[ir_now] += jB_divided_IB * dlambda_dne * ((-1.0 / (lambda_s * lambda_s)) * survivalProb + 1 / lambda_s * I_s_squared * survivalProb) * dV_beam_prime;
            }
        }
    }
}

/**
 * Set the Jacobian elements corresponding to this term.
 */
bool NBIElectronHeatTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId,
    FVM::Matrix *jac, const real_t *x){
    if (derivId != id_ncold && derivId != id_Tcold &&
        derivId != id_ion_density && derivId != id_ion_temperature)
        return false; 

    for (len_t ir = 0; ir < nr; ++ir){
        real_t dP = 0.0;

        if (derivId == id_ncold)
            dP = H_r_dne[ir];
        else if (derivId == id_Tcold)
            dP = H_r_dTe[ir];
        else if (derivId == id_ion_density)
            dP = H_r_dni[ir];
        else if (derivId == id_ion_temperature)
            dP = H_r_dTi[ir];

        jac->SetElement(ir, ir, dP);
    }

    return false;
}

/**
 * Set the non-linear function vector for this term.
 */
void NBIElectronHeatTerm::SetVectorElements(real_t *rhs, const real_t *x){
    for (len_t ir = 0; ir < nr; ir++){
        rhs[ir] -= NBIHeatTerm[ir];
    }
}

/**
 * Sets the matrix elements for this term.
 */
void NBIElectronHeatTerm::SetMatrixElements(FVM::Matrix *mat, real_t *rhs){
    const real_t *n_cold = unknowns->GetUnknownData(id_ncold);
    for (len_t ir = 0; ir < nr; ++ir){
        real_t factor = 2.0 / (3.0 * n_cold[ir]);
        rhs[ir] -= factor * NBIHeatTerm[ir];
    }
}
