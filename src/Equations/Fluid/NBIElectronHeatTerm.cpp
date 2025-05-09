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

using namespace DREAM;
/**
 * Constructor
 */
NBIElectronHeatTerm::NBIElectronHeatTerm(FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns, ADAS *adas,
                                         real_t s_max, real_t r_beam,
                                         const real_t P0[3], const real_t n[3],
                                         real_t Ti_beam, real_t m_i_beam,
                                         real_t beamPower,
                                         FVM::Interpolator1D *j_B_profile, real_t Z0, real_t Zion, real_t R0) : EquationTerm(grid)
{

    // Set grid
    this->radialGrid = grid->GetRadialGrid();
    this->nr = grid->GetNr();
    this->NBIHeatTerm = new real_t[nr];
    this->Deposition_profile = new real_t[nr];
    this->Deposition_profile_times_Vprime = new real_t[nr];

    // Get from other classes
    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_ion_density = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    this->id_ion_temperature = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->adas = adas;

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
    this->beamPower = beamPower;
    this->Z0 = Z0;
    this->Zion = Zion;
    this->R0 = R0;
    for (int i = 0; i < 3; ++i)
    {
        this->P0[i] = P0[i];
        this->n[i] = n[i];
    }
    this->s_start = 0;
    this->j_B_profile = j_B_profile;
    this->n_beam_radius = 50;
    this->d_beam_radius = r_beam / n_beam_radius;
    this->n_beam_theta = 50;
    this->d_beam_theta = 2.0 * M_PI / n_beam_theta;
    this->n_beam_s = 50;
    this->d_beam_s = s_max / n_beam_s;
    PrecomputeBeamBasisVectors();
    PrecomputeFluxSurfaces();

    // Integrate VpVol to get plasma volume
    this->plasmaVolume = 0;
    const real_t *VpVol = grid->GetVpVol();
    const real_t *dr = this->radialGrid->GetDr();
    for (len_t ir = 0; ir < nr; ir++)
        this->plasmaVolume += VpVol[ir] * dr[ir] * R0;

    // Trapezoidal integration to get total beam current
    this->I_B = 0.0;
    const real_t *r_beam_grid = this->j_B_profile->GetX();
    len_t N = this->j_B_profile->GetNx();
    for (len_t i = 0; i < N - 1; ++i)
    {
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
 * Precomputes the flux surfaces in R and Z.
 */
void NBIElectronHeatTerm::PrecomputeFluxSurfaces()
{
    cachedFluxSurfaces.resize(nr);
    for (len_t ir = 0; ir < nr; ++ir)
    {
        cachedFluxSurfaces[ir].resize(ntheta);
        for (len_t itheta = 0; itheta < ntheta; ++itheta)
        {
            real_t theta = itheta * dtheta;
            real_t R = R0 + radialGrid->GetFluxSurfaceRMinusR0_theta(ir, theta);
            real_t Z = radialGrid->GetFluxSurfaceZMinusZ0_theta(ir, theta);
            cachedFluxSurfaces[ir][itheta] = {R, Z};
        }
    }
}

/**
 * Precomputes the beam basis vectors of beam e1 and e2.
 */
void NBIElectronHeatTerm::PrecomputeBeamBasisVectors()
{
    real_t norm_n = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (int i = 0; i < 3; ++i)
        n_norm[i] = n[i] / norm_n;

    std::array<real_t, 3> a;
    if (std::abs(n_norm[2]) < 0.9)
    {
        a = {0, 0, 1};
    }
    else
    {
        a = {1, 0, 0};
    }
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
NBIElectronHeatTerm::~NBIElectronHeatTerm()
{
    delete[] NBIHeatTerm;
    delete[] Deposition_profile;
    delete[] Deposition_profile_times_Vprime;
}

/**
 * Rebuild: Called at the start of every timestep.
 * Compute added heating term for each timestep and flux surface radius
 */
void NBIElectronHeatTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *unknowns)
{
    for (len_t ir = 0; ir < nr; ++ir)
    {
        NBIHeatTerm[ir] = 0.0;
        Deposition_profile[ir] = 0.0;
        Deposition_profile_times_Vprime[ir] = 0.0;
    }
    // Main Calculation
    ComputeDepositionProfile(unknowns);
    // Compute the deposited fraction of the beam
    real_t cumulative = 0.0;
    real_t frac = 0.0;
    for (len_t ir = 0; ir < nr; ir++)
    {
        cumulative += Deposition_profile_times_Vprime[ir] * d_beam_radius;
        frac = cumulative;
    }
    printf("✅ Deposited fraction of beam = %.6f\n", cumulative);
    // Compute volume-weighted integral of raw deposition profile
    real_t integral = 0.0;
    for (len_t ir = 0; ir < nr; ++ir)
        integral += Deposition_profile[ir] * this->radialGrid->GetVpVol()[ir] * R0 * this->radialGrid->GetDr()[ir];

    // Now normalize so that ∫ H(r) V'(r) dr = plasmaVolume
    if (integral > 0)
    {
        real_t norm_factor = plasmaVolume / integral;
        for (len_t ir = 0; ir < nr; ++ir)
            Deposition_profile[ir] *= norm_factor;
    }

    for (len_t ir = 0; ir < nr; ir++)
    {
        NBIHeatTerm[ir] = Deposition_profile[ir] * beamPower / plasmaVolume * frac;
    }
}

/**
 * Computes and mean free path variable for the NBI heating term
 */
real_t NBIElectronHeatTerm::ComputeMeanFreePath(len_t ir)
{
    // Ionization rate coefficient (SCD)
    real_t ncold = this->unknowns->GetUnknownData(id_ncold)[ir];
    real_t Tcold = this->unknowns->GetUnknownData(id_Tcold)[ir];
    real_t ni = this->unknowns->GetUnknownData(id_ion_density)[ir];
    real_t Ti = this->unknowns->GetUnknownData(id_ion_temperature)[ir];

    ADASRateInterpolator *scd = adas->GetSCD(Zion);
    ADASRateInterpolator *ccd = adas->GetCCD(Zion);

    real_t I_ion = scd->Eval(Z0, ncold, Tcold);
    real_t I_CX = ccd->Eval(Z0, ni, Ti);

    real_t v_NBI = std::sqrt(2.0 * Ti_beam / m_i_beam);
    real_t MeanFreePath = v_NBI / (ncold * (I_ion + I_CX));

    return MeanFreePath;
}

/**
 * Computes the flux surface index for a given pencil beam position (s_B, r_B, theta_B).
 * This function checks if the pencil beam intersects with any flux surface and returns the index of the first one it finds.
 */
int_t NBIElectronHeatTerm::CalculatePencilBeamFindFlux(real_t s_B, real_t r_B, real_t theta_B)
{
    // Convert the beam coordinates (s_B, r_B, theta_B) to Cartesian coordinates
    real_t dx = r_B * cos(theta_B);
    real_t dy = r_B * sin(theta_B);
    std::array<real_t, 3> P0_current = {
        P0[0] + dx * e1[0] + dy * e2[0],
        P0[1] + dx * e1[1] + dy * e2[1],
        P0[2] + dx * e1[2] + dy * e2[2]};
    int_t best_ir = -1;
    real_t smallest_diff = std::numeric_limits<real_t>::infinity();

    for (len_t ir = 0; ir < nr; ir++)
    { // Loop through all flux surfaces
        for (len_t i = 0; i < ntheta; i++)
        { // Loop thourgh theta and discretizise
            len_t j = (i + 1) % ntheta;

            const FluxSurfacePoint &pt1 = cachedFluxSurfaces[ir][i];
            const FluxSurfacePoint &pt2 = cachedFluxSurfaces[ir][j];

            real_t R1 = pt1.R, Z1 = pt1.Z;
            real_t R2 = pt2.R, Z2 = pt2.Z;

            real_t dZ = Z2 - Z1;
            real_t dR = R2 - R1;

            real_t ratio = dR / dZ;

            real_t a0 = R1 * R1 - P0_current[0] * P0_current[0] - P0_current[1] * P0_current[1] + 2 * R1 * (P0_current[2] - Z1) * ratio + (P0_current[2] - Z1) * (P0_current[2] - Z1) * ratio * ratio;

            real_t a1 = 2 * R1 * n[2] * ratio + 2 * n[2] * (P0_current[2] - Z1) * ratio * ratio - 2 * (P0_current[0] * n[0] + P0_current[1] * n[1]);

            real_t a2 = n[2] * n[2] * ratio * ratio - n[0] * n[0] - n[1] * n[1];

            real_t discriminant = a1 * a1 / (4 * a2 * a2) - a0 / a2;
            if (discriminant < 0)
                continue;

            real_t sqrt_disc = std::sqrt(discriminant);
            real_t l1 = (-a1 / (2 * a2) + sqrt_disc);
            real_t l2 = (-a1 / (2 * a2) - sqrt_disc);

            // Check both solutions
            auto check_solution = [&](real_t l) -> std::pair<bool, real_t>
            {
                if (l < 0)
                    return {false, 0};
                real_t t = (P0_current[2] - Z1 + l * n[2]) / dZ;
                if (t < 0 || t > 1)
                    return {false, 0};
                return {true, std::abs(l - s_B)};
            };

            // Check both solutions and update if better match found
            auto [valid1, diff1] = check_solution(l1);
            if (valid1 && diff1 < smallest_diff)
            {
                smallest_diff = diff1;
                best_ir = ir;
            }

            auto [valid2, diff2] = check_solution(l2);
            if (valid2 && diff2 < smallest_diff)
            {
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
    const real_t a[3],
    real_t &r, real_t &theta, real_t &s)
{

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
void NBIElectronHeatTerm::ComputeDepositionProfile(
    FVM::UnknownQuantityHandler *unknowns)
{
    for (len_t ir = 0; ir < nr; ++ir)
    {
        Deposition_profile_times_Vprime[ir] = 0.0;
        Deposition_profile[ir] = 0.0;
    }

    for (real_t i_beam_theta = 0; i_beam_theta < n_beam_theta; i_beam_theta++)
    {
        real_t beam_theta = i_beam_theta * d_beam_theta;
        for (real_t i_beam_radius = 0; i_beam_radius < n_beam_radius; i_beam_radius++)
        {
            real_t beam_radius = i_beam_radius * d_beam_radius;
            real_t I_s = 0.0;
            for (real_t i_beam_s = 0; i_beam_s < n_beam_s; i_beam_s++)
            {
                real_t beam_s = s_start + i_beam_s * d_beam_s;

                len_t ir_now = CalculatePencilBeamFindFlux(beam_s, beam_radius, beam_theta);
                if (ir_now < 0 || ir_now >= nr)
                {
                    printf("⚠️  Warning: ir_now = %d out of bounds (s_B = %.6e, r_B = %.6e, theta_B = %.6e)\n", ir_now, beam_s, beam_radius, beam_theta);
                    continue;
                }
                real_t lambda_s = ComputeMeanFreePath(ir_now);

                I_s += 1 / lambda_s * d_beam_s;
                real_t survivalProb = exp(-I_s);

                const real_t *j_B_values = j_B_profile->Eval(beam_radius);
                real_t j_B_r = j_B_values[0];

                real_t dV_beam_prime = beam_radius * d_beam_theta * d_beam_s;
                real_t H_r_theta_phi = (j_B_r / I_B) * (1.0 / lambda_s) * survivalProb;
                Deposition_profile[ir_now] += H_r_theta_phi;
                Deposition_profile_times_Vprime[ir_now] += H_r_theta_phi * dV_beam_prime;
            }
        }
    }
}

/**
 * Set the Jacobian elements corresponding to this term.
 */
bool NBIElectronHeatTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId,
    FVM::Matrix *jac, const real_t *x)
{
    for (len_t ir = 0; ir < nr; ir++)
    {
        real_t dP = Compute_dP_derivative(ir, derivId, unknowns);
        jac->SetElement(ir, ir, dP);
    }
    return true;
}

/**
 * Set the non-linear function vector for this term.
 */
void NBIElectronHeatTerm::SetVectorElements(real_t *rhs, const real_t *x)
{
    for (len_t ir = 0; ir < nr; ir++)
    {
        rhs[ir] += NBIHeatTerm[ir];
    }
}

/**
 * Sets the matrix elements for this term.
 */
void NBIElectronHeatTerm::SetMatrixElements(FVM::Matrix *mat, real_t *rhs)
{
    if (rhs == nullptr)
        return;
    const real_t *n_cold = unknowns->GetUnknownData(id_ncold);
    for (len_t ir = 0; ir < nr; ir++)
    {
        real_t factor = 2.0 / (3.0 * n_cold[ir]);
        rhs[ir] += factor * NBIHeatTerm[ir];
    }
}
