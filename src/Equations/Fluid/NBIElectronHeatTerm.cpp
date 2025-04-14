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


namespace DREAM{
/**
 * Constructor
 */
NBIElectronHeatTerm::NBIElectronHeatTerm(FVM::Grid* grid, FVM::UnknownQuantityHandler *unknowns, ADAS* adas, 
    real_t ds, real_t s_max, real_t r_beam,
    std::array<real_t, 3> P0, std::array<real_t, 3> n,
    real_t Ti_beam, real_t m_i_beam,
    real_t beamPower, real_t plasmaVolume,
    FVM::Interpolator1D *j_B_profile, real_t Z0, real_t Zion, real_t R0) : EquationTerm(grid) {
    
    //Set grid    
    this->radialGrid = grid->GetRadialGrid(); 
    this->nr = grid->GetNr();
    this->NBIHeatTerm = new real_t[nr];

    //Get from other classes
    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_ion_density = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS); //UQTY_ION_SPECIES
    this->id_ion_temperature = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER); 
    this->adas = adas;

    // Tell DREAM these are dependencies
    this->AddUnknownForJacobian(unknowns, id_ncold);
    this->AddUnknownForJacobian(unknowns, id_Tcold);
    this->AddUnknownForJacobian(unknowns, id_ion_density);
    this->AddUnknownForJacobian(unknowns, id_ion_temperature);
    
    //Set integration step size
    this->ds = ds;
    this->unknowns = unknowns;
    this-> ntheta = 100;
    this-> dtheta = 2.0 * M_PI / ntheta;
    this-> nphi = 100;
    this-> dphi = 2.0*M_PI/nphi;
    
    //Set beam parameters
    this->s_max = s_max;
    this->r_beam = r_beam;
    this->P0 = P0;
    this->n = n;
    this->Ti_beam = Ti_beam;
    this->m_i_beam = m_i_beam;
    this->beamPower = beamPower;
    this->plasmaVolume = plasmaVolume;
    this->Z0 = Z0;    
    this->Zion = Zion;  
    this->R0 = R0;

    this->j_B_profile = j_B_profile;
    //Precompute the flux surfaces and the beam basis vectors
    PrecomputeFluxSurfaces();
    PrecomputeBeamBasisVectors(); 
    
    
    printf("NBIElectronHeatTerm: Initialized with parameters:\n");
    printf("  ds = %e\n", ds);
    printf("  s_max = %e\n", s_max);
    printf("  r_beam = %e\n", r_beam);  

    // Calculate I_B as integral of j_B over radius
    this->I_B = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const real_t* j_B_values = j_B_profile->Eval(ir);
        real_t dr = (ir < nr-1) ? //diff calc
            radialGrid->GetR(ir+1) - radialGrid->GetR(ir) : 
            radialGrid->GetR(ir) - radialGrid->GetR(ir-1);
        this->I_B += j_B_values[0] * dr;  // Integrate j_B over radius
    }

    printf("  Calculated I_B = %e\n", this->I_B);
}


/**
 * Precomputes the flux surfaces in R and Z.
 */
void NBIElectronHeatTerm::PrecomputeFluxSurfaces() {
    cachedFluxSurfaces.resize(nr);
    for (len_t ir = 0; ir < nr; ir++) {
        cachedFluxSurfaces[ir].resize(ntheta);
        for (len_t itheta = 0; itheta < ntheta; itheta++) {
            real_t theta = itheta * dtheta;
            real_t R = R0 + radialGrid->GetFluxSurfaceR(ir, theta); 
            real_t Z = radialGrid->GetFluxSurfaceZ(ir, theta);
            cachedFluxSurfaces[ir][itheta] = { R, Z };
        }
    }
}

/**
 * Precomputes the beam basis vectors e1 and e2.
 */
void NBIElectronHeatTerm::PrecomputeBeamBasisVectors() {
    // Normalize the beam direction vector 'n'
    real_t norm_n = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    for (int i = 0; i < 3; ++i)
        n_norm[i] = n[i] / norm_n;

    // Choose 'a' to be perpendicular to the beam direction
    std::array<real_t, 3> a;
    if (std::abs(n_norm[2]) < 0.9) {
        // If beam is not mostly vertical, use vertical direction
        a = {0, 0, 1};
    } else {
        // If beam is mostly vertical, use horizontal direction
        a = {1, 0, 0};
    }

    // Compute e1 (perpendicular to beam direction)
    real_t a_dot_n = a[0]*n_norm[0] + a[1]*n_norm[1] + a[2]*n_norm[2];
    for (int i = 0; i < 3; ++i)
        e1[i] = a[i] - a_dot_n * n_norm[i];

    real_t norm_e1 = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
    for (int i = 0; i < 3; ++i)
        e1[i] /= norm_e1;
    
    // Compute e2 as the cross product of n_norm and e1
    e2[0] = n_norm[1]*e1[2] - n_norm[2]*e1[1];
    e2[1] = n_norm[2]*e1[0] - n_norm[0]*e1[2];
    e2[2] = n_norm[0]*e1[1] - n_norm[1]*e1[0];
}

/**
 * Destructor
 */
NBIElectronHeatTerm::~NBIElectronHeatTerm() {
    delete[] NBIHeatTerm;
}

/**
 * Rebuild: Called at the start of every timestep.
 * Compute added heating term for each timestep and flux surface radius
 */
void NBIElectronHeatTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler* unknowns)
{
    printf("NBIElectronHeatTerm::Rebuild called.\n");
    for (len_t ir = 0; ir < nr;ir ++){
        real_t ncold = unknowns ->GetUnknownData(id_ncold)[ir];
        real_t Tcold = unknowns ->GetUnknownData(id_Tcold)[ir];
        real_t ni = unknowns ->GetUnknownData(id_ion_density)[ir];
        real_t Ti = unknowns ->GetUnknownData(id_ion_temperature)[ir];

        printf("Checking flux surface %llu...\n", ir);
        NBIHeatTerm[ir] = beamPower/plasmaVolume *ComputeDepositionProfile(ir,ncold,Tcold,Ti,ni, unknowns); //compute the NBI heating term * times the current density and beam intensisty
        printf("[ir=%llu] Deposition = %e\n", ir, NBIHeatTerm[ir]);

    }
}


/**
 * Computes and intigrates the Mean Free Path variable for the NBI heating term
 */
real_t NBIElectronHeatTerm::ComputeMeanFreePath(len_t ir, real_t ncold, real_t Tcold, real_t Ti, real_t ni) {
    // Ionization rate coefficient (SCD)
    ADASRateInterpolator* scd = adas->GetSCD(Zion);
    // Charge-exchange rate coefficient (CCD)
    ADASRateInterpolator* ccd = adas->GetCCD(Zion);
    //Evaluate
    real_t I_ion = scd->Eval(Z0, ncold, Tcold);
    real_t I_CX = ccd->Eval(Z0, ni, Ti);
    real_t v_NBI = sqrt(2.0 * Ti_beam / m_i_beam); 

    real_t MeanFreePath = v_NBI /(ncold * (I_ion+I_CX));
    return MeanFreePath; 
}

/**
 * Computes the flux surface index for a given pencil beam position (s_B, r_B, theta_B).
 * This function checks if the pencil beam intersects with any flux surface and returns the index of the first one it finds.
 */
len_t NBIElectronHeatTerm::CalculatePencilBeamFindFlux(real_t s_B, real_t r_B, real_t theta_B) {
    // Convert the beam coordinates (s_B, r_B, theta_B) to Cartesian coordinates
    real_t dx = r_B * cos(theta_B);
    real_t dy = r_B * sin(theta_B);
    std::array<real_t, 3> P0_current = {
        P0[0] + dx*e1[0] + dy*e2[0],
        P0[1] + dx*e1[1] + dy*e2[1],
        P0[2] + dx*e1[2] + dy*e2[2]
    };

    const real_t tol_s = 2e-1;
    for (len_t ir = 0; ir < nr; ir++) { //Loop through all flux surfaces
        for (len_t i = 0; i < ntheta; i++) { //Loop thourgh theta and discretizise
            len_t j = (i + 1) % ntheta;

            const FluxSurfacePoint& pt1 = cachedFluxSurfaces[ir][i];
            const FluxSurfacePoint& pt2 = cachedFluxSurfaces[ir][j];

            real_t R1 = pt1.R, Z1 = pt1.Z;
            real_t R2 = pt2.R, Z2 = pt2.Z;

            real_t dZ = Z2 - Z1;
            real_t dR = R2 - R1;

            if (std::abs(dZ) < 1e-12)
                continue;

            real_t ratio = dR / dZ;

            real_t a0 = R1*R1 - P0_current[0]*P0_current[0] - P0_current[1]*P0_current[1]
                      + 2*R1*(P0_current[2]-Z1)*ratio
                      + (P0_current[2]-Z1)*(P0_current[2]-Z1)*ratio*ratio;

            real_t a1 = 2*R1*n[2]*ratio
                      + 2*n[2]*(P0_current[2]-Z1)*ratio*ratio
                      - 2*(P0_current[0]*n[0] + P0_current[1]*n[1]);

            real_t a2 = n[2]*n[2]*ratio*ratio - n[0]*n[0] - n[1]*n[1];

            if (std::abs(a2) < 1e-12)
                continue;

            real_t discriminant = a1*a1/(4*a2*a2) - a0/a2;
            if (discriminant < 0)
                continue;

            real_t sqrt_disc = sqrt(discriminant);
            real_t l1 = (-a1/(2*a2) + sqrt_disc);
            real_t l2 = (-a1/(2*a2) - sqrt_disc);
            
            // Help function to check if the solutions are valid
            auto check_solution = [&](real_t l) -> bool {
                if (l < 0) return false;
                real_t t = (P0_current[2] - Z1 + l * n[2]) / dZ;
                return (t >= 0 && t <= 1);
            };

            // Check if the solutions are valid
            std::vector<real_t> valid_s;
            if (check_solution(l1)) valid_s.push_back(l1);
            if (check_solution(l2)) valid_s.push_back(l2);
            
            // Check if any valid solution is close to s_B
            for (real_t s_val : valid_s) {
                if (std::fabs(s_val - s_B) < tol_s) {
                    return ir;
                }
                
            }
        }
    }
    
    return static_cast<len_t>(-1);
}



/**
 * Computes the survival probabiltiy at one flux surface point by intigrating the correspodning pencil beam up to that point 
 */
real_t NBIElectronHeatTerm::ComputeSurvivalProbability(real_t s_B, real_t r_B, real_t theta_B, FVM::UnknownQuantityHandler* unknowns) 
{
    std::vector<real_t> lambda_cache(nr, -1.0);  // -1 = not computed yet
    real_t I_s = 0.0; 
    real_t s_step = 0.0;

    while (s_step <= s_B) {
        len_t ir_now = CalculatePencilBeamFindFlux(s_step, r_B, theta_B);
        if (ir_now == static_cast<len_t>(-1)) {
        s_step += ds;  
        continue;}

        real_t ncold_now = unknowns->GetUnknownData(id_ncold)[ir_now];
        real_t Tcold_now = unknowns->GetUnknownData(id_Tcold)[ir_now];
        real_t Ti_now = unknowns->GetUnknownData(id_ion_temperature)[ir_now];
        real_t ni_now = unknowns->GetUnknownData(id_ion_density)[ir_now];
        
        //Check if the values has been calculated and stored before
        real_t lambda_now;
        if (lambda_cache[ir_now] < 0) {
            lambda_now = ComputeMeanFreePath(ir_now, ncold_now, Tcold_now, Ti_now, ni_now);
            lambda_cache[ir_now] = lambda_now;
        } else {
            lambda_now = lambda_cache[ir_now];
        }

        I_s += (ds / lambda_now);
        s_step += ds;  
    }
    if (I_s == 0) {
    printf("⚠️  Warning: I_s = 0 → No contribution to survival probability (s_B = %.6e), (r_B = %.6e), (theta_B = %.6e). Some part of the beam is probabily outside the plasma \n", s_B, r_B, theta_B);
    return 0;  
}
    return exp(-I_s);  // Survival probability
}

/**
 * Takes in the beam paramaters and turns a point on a flux surface (in x y z) into the beam coordinate system 
 */
void NBIElectronHeatTerm::CartesianToCylindrical(
    real_t x, real_t y, real_t z,
    const std::array<real_t, 3>& P0, 
    const std::array<real_t, 3>& n,  
    const std::array<real_t, 3>& a,  
    real_t &r, real_t &theta, real_t &s
) {
    
   //Convert cartesian coordinates to beamline coordinates
    std::array<real_t, 3> d = {x - P0[0], y - P0[1], z - P0[2]};
    s = d[0]*n_norm[0] + d[1]*n_norm[1] + d[2]*n_norm[2]; //s

    std::array<real_t, 3> d_trans;
    for (int i = 0; i < 3; ++i)
        d_trans[i] = d[i] - s * n_norm[i];

    real_t r_cos_theta = d_trans[0]*e1[0] + d_trans[1]*e1[1] + d_trans[2]*e1[2];
    real_t r_sin_theta = d_trans[0]*e2[0] + d_trans[1]*e2[1] + d_trans[2]*e2[2];
    r = sqrt(r_cos_theta*r_cos_theta + r_sin_theta*r_sin_theta);
    theta = atan2(r_sin_theta, r_cos_theta); //theta
    if (theta < 0)
        theta += 2 * M_PI; //Get the corresponding positive angle
}

/**
 * Computes the deposition profile for a flux surface
 */
real_t NBIElectronHeatTerm::ComputeDepositionProfile(len_t ir, real_t ncold, real_t Tcold, real_t Ti, real_t ni, FVM::UnknownQuantityHandler* unknowns
){
    real_t total_deposition = 0.0;
    real_t volume_element = 0.0;
    real_t lambda_s = ComputeMeanFreePath(ir,ncold,Tcold,Ti,ni); // At current flux surface point
    
    for (len_t itheta = 0; itheta < ntheta; itheta++){ //Integrate by theta
        const FluxSurfacePoint& pt = cachedFluxSurfaces[ir][itheta];
        real_t R = pt.R;
        real_t Z = pt.Z;
        //Compute Configuration space Jacobian
        real_t theta = itheta*dtheta;
        
        for (len_t iphi = 0; iphi < nphi; iphi++){ //Integrate by phi
            //For every point on that flux surface
            real_t phi = iphi*dphi;
            real_t x = R * cos(phi); 
            real_t y = R * sin(phi);
            real_t z = Z;
        
            // Compute the beamline coordinates r_B,theta_B,s_B for this/every point, so now we are in the beam system
            real_t r_B, theta_B, s_B;
            CartesianToCylindrical(x, y, z, P0, n, a, r_B, theta_B, s_B);

            //Check if the point is within the beam radius and beam length
            if (r_B >= r_beam) 
                continue;
            if (s_B < 0 || s_B > s_max)
                continue;
            //Compute survival probability at this point
            real_t survivalProb = ComputeSurvivalProbability(s_B,r_B,theta_B, unknowns);
            const real_t* j_B_values = j_B_profile->Eval(ir);
            real_t j_B_r = j_B_values[0];
            real_t J = this->radialGrid->ComputeConfigurationSpaceJacobian(ir, theta); 
            // Integrate the deposition profile over the flux surface
            real_t H_r_theta_phi = (j_B_r / I_B) * (1.0 / lambda_s) * survivalProb;
            total_deposition += H_r_theta_phi * J * dtheta * dphi; 
            volume_element += J * dtheta * dphi; 
        }
    }
    // If the deposition or volume element is NaN or zero, return 0
    if (std::isnan(total_deposition) || std::isnan(volume_element) || volume_element == 0) {
        printf("⚠️  Warning: Invalid deposition result for flux surface %lu (total_deposition = %.3f, volume_element = %.3f). The beam is not depositing at this flux surface.\n", ir, total_deposition, volume_element);
        return 0.0;
    }
    return total_deposition / volume_element;
}

/**
 * Set the Jacobian elements corresponding to this term.
 */
bool NBIElectronHeatTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId,
    FVM::Matrix *jac, const real_t *x
)
 {
    for (len_t ir = 0; ir< nr; ir++) {
        auto [dP_dni, dP_dTi, dP_dTe, dP_dne] = Compute_dP_derivative(ir, unknowns);

        real_t dP_dX = 0.0;
        if (derivId == id_ion_density)
            dP_dX = dP_dni;
        else if (derivId == id_ion_temperature)
            dP_dX = dP_dTi;
        else if (derivId == id_ncold)
            dP_dX = dP_dne;
        else if (derivId == id_Tcold)
            dP_dX = dP_dTe;

        // Apply to Jacobian matrix (NBI heating appears in the equation for T_cold)
        jac->SetElement(id_Tcold, ir, dP_dX);

    }

    return true;
}



/**
 * Set the non-linear function vector for this term.
 */
void NBIElectronHeatTerm::SetVectorElements(real_t *rhs, const real_t *x) {
    for (len_t ir = 0; ir < nr; ir++){
        rhs[ir] += NBIHeatTerm[ir];
        printf("NBIHeatTerm[%d] = %e\n", ir, NBIHeatTerm[ir]);
        }

}

/**
 * Sets the matrix elements for this term.
 */

 void NBIElectronHeatTerm::SetMatrixElements(FVM::Matrix *mat, real_t *rhs) {
    if (rhs == nullptr)
        return;
    const real_t *n_cold = unknowns->GetUnknownData(id_ncold);
    
    // Add heating contribution to each radial point
    for (len_t ir = 0; ir < nr; ir++) {
        real_t factor = 2.0 / (3.0 * n_cold[ir]);
        rhs[ir] += factor * NBIHeatTerm[ir];
    }
}


}


















