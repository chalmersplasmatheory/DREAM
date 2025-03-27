/**
 * Implementation of equation term to the heating of cold electrons from an NBI beam 
 */
#include <cmath>
#include <array>
#include "DREAM/Equations/Fluid/NBIElectronHeatTerm.hpp" 
#include "FVM/Grid/RadialGrid.hpp"
#include "DREAM/ADAS.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Equations/Fluid/IonRateEquation.hpp"
#include <tuple>
#include <algorithm>
#include "FVM/UnknownQuantityHandler.hpp"






namespace DREAM{
/**
 * Constructor
 */
NBIElectronHeatTerm::NBIElectronHeatTerm(FVM::Grid* grid, FVM::UnknownQuantityHandler *unknowns, ADAS* adas, 
    real_t ds, real_t s_max, real_t r_beam,
    std::array<real_t, 3> P0, std::array<real_t, 3> n,
    real_t Ti_beam, real_t m_i_beam,
    real_t beamPower, real_t plasmaVolume,
    real_t j_B, real_t I_B) : EquationTerm(grid) {

    this->NBIHeatTerm = new real_t[nr];
    //Get from other classes
    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_ion_density = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this->id_ion_temperature = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER); 

    // Tell DREAM these are dependencies
    this->AddUnknownForJacobian(unknowns, id_ncold);
    this->AddUnknownForJacobian(unknowns, id_Tcold);
    this->AddUnknownForJacobian(unknowns, id_ion_density);
    this->AddUnknownForJacobian(unknowns, id_ion_temperature);



    this->radialGrid = grid->GetRadialGrid(); //we have not defined whcih grid to use
    this->nr = grid->GetNr();

    this->ds = ds;
    this->s_max = s_max;
    this->r_beam = r_beam;
    this->P0 = P0;
    this->n = n;
    this->Ti_beam = Ti_beam;
    this->m_i_beam = m_i_beam;
    this->beamPower = beamPower;
    this->plasmaVolume = plasmaVolume;
    this->j_B = j_B;
    this->I_B = I_B;

    //Check that not beam direction so GS-work
    this-> a  = {0, 0, 1}; 
    if (n == a)
        this-> a = {0, 1, 0};

    //For the mean free path
    this-> Z0 = 0;    // ex
    this-> Zion = 1;  // ex

    PrecomputeBeamIntersections(P0,n); 
    
}
/**
 * Pre-Computes all the intersections between the beam and the fluxsurface. 
 This to be used when intigrating along the beam and mapping each step to a step on the beam.
 */

void NBIElectronHeatTerm::PrecomputeBeamIntersections(const std::array<real_t, 3>& P0,
    std::array<real_t, 3>& n) {
    fluxSurfaceIntersectionPoints.clear(); //All flux surcafe intersections for total beam
    fluxSurfaceToS.clear(); //Maps intersections to a flux surface

    real_t norm_n = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    for (int i = 0; i < 3; ++i)
        n[i] /= norm_n;

    const len_t ntheta = 200; 
    const real_t dtheta = 2 * M_PI / ntheta;
    //Walks though every flux surface and find the intersection point(s).
    for (len_t ir = 0; ir < this->nr; ir++) {
        for (len_t it = 0; it < ntheta - 1; it++) {
            real_t theta1 = it * dtheta;
            real_t theta2 = (it + 1) * dtheta;

            real_t x1 = radialGrid->GetFluxSurfaceR(ir, theta1);
            real_t z1 = radialGrid->GetFluxSurfaceZ(ir, theta1);
            real_t x2 = radialGrid->GetFluxSurfaceR(ir, theta2);
            real_t z2 = radialGrid->GetFluxSurfaceZ(ir, theta2);

            real_t dz = z2 - z1;
            real_t dxi = x2 - x1;
            real_t zi_z0 = z1 - P0[2];
            if (std::abs(dz) < 1e-12)
                continue; //avoid 0 
            real_t frac = dxi / dz;
            real_t a2 = n[2]*n[2]*frac*frac - (n[0]*n[0] + n[1]*n[1]);
            real_t a1 = 2 * n[0]*n[2]*frac;
            real_t a0 = x1*x1 - P0[0]*P0[0] - P0[1]*P0[1] +
                        2 * x1 * zi_z0 * frac +
                        2 * n[2] * zi_z0 * frac +
                        zi_z0*zi_z0 * frac*frac;
            real_t discriminant = a1*a1 - 4*a2*a0;
            if (discriminant < 0)
                continue; //avoid 0 

            real_t sqrt_disc = sqrt(discriminant);
            for (int sign = -1; sign <= 1; sign += 2) 
            {
                real_t l = (-a1 + sign * sqrt_disc) / (2 * a2);
                if (l < 0)
                    continue;//avoid 0 
                
                real_t z_b = P0[2] + l * n[2];
                real_t t = (z_b - z1) / dz;
                if (t < 0 || t > 1)
                    continue;//check if satisfies t

                fluxSurfaceIntersectionPoints.push_back(l); 
                fluxSurfaceToS[ir].push_back(l); 
            }
        }
    } 

    std::sort(fluxSurfaceIntersectionPoints.begin(), fluxSurfaceIntersectionPoints.end()); //Sort in length-order

    // Sort the values for each key (flux surface ir)
    for (auto& kv : fluxSurfaceToS)
        std::sort(kv.second.begin(), kv.second.end());
}

/**
 * Help function that returns the flux surface index ir for a given beamline position s.
 * This performs a lookup in the precomputed fluxSurfaceToS map.
 */
len_t NBIElectronHeatTerm::FindFluxSurfaceIndexForS(real_t s) {
    for (const auto& [ir, s_list] : fluxSurfaceToS) { //Get Key,Value from fluxSurfaceToS
        for (size_t i = 0; i + 1 < s_list.size(); i += 2) { //Go though pair
            if (s >= s_list[i] && s <= s_list[i + 1]) //Check if in pair
                return ir;
        }
    }
    return static_cast<len_t>(-1);  // No find
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

    for (len_t ir = 0; ir < nr;ir ++){
        real_t ncold = unknowns ->GetUnknownData(id_ncold)[ir];
        real_t Tcold = unknowns ->GetUnknownData(id_Tcold)[ir];
        real_t ni = unknowns ->GetUnknownData(id_ion_density)[ir];
        real_t Ti = unknowns ->GetUnknownData(id_ion_temperature)[ir];

        
        NBIHeatTerm[ir] = beamPower/plasmaVolume *ComputeDepositionProfile(ir,ncold,Tcold,Ti,ni, unknowns); //compute the NBI heating term * times the current density and beam intensisty
    }
}


/**
 * Computes and intigrates  Mean Free Path
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
 * Computes the survival probabiltiy at one flux surface point by intigrating the correspodning pencil beam up to that point 
 */
real_t NBIElectronHeatTerm::ComputeSurvivalProbability(real_t ir, real_t s_B, FVM::UnknownQuantityHandler* unknowns) {
    //For a flux surface point, but expressed in cylindrical/beam coordinates
    real_t I_s = 0.0; 
    real_t s = 0.0;
    while (s < s_B) {
        // Find the corresponding flux surface for the current beam position
        len_t ir_now = FindFluxSurfaceIndexForS(s);
        if (ir_now == static_cast<len_t>(-1)) //Take previos one if do not find
            ir_now = nr - 1;

        real_t ncold_now = unknowns->GetUnknownData(id_ncold)[ir_now];
        real_t Tcold_now = unknowns->GetUnknownData(id_Tcold)[ir_now];
        real_t Ti_now = unknowns->GetUnknownData(id_ion_temperature)[ir_now];
        real_t ni_now = unknowns->GetUnknownData(id_ion_density)[ir_now];


        real_t lambda_now = ComputeMeanFreePath(ir_now, ncold_now, Tcold_now, Ti_now, ni_now);
        I_s += (ds / lambda_now);
        s += ds;  
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
    std::array<real_t, 3> n_norm = n;
    real_t norm_n = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    for (int i = 0; i < 3; ++i)
        n_norm[i] /= norm_n;

    std::array<real_t, 3> e1;
    real_t a_dot_n = a[0]*n_norm[0] + a[1]*n_norm[1] + a[2]*n_norm[2];
    for (int i = 0; i < 3; ++i)
        e1[i] = a[i] - a_dot_n * n_norm[i];

    real_t norm_e1 = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
    for (int i = 0; i < 3; ++i)
        e1[i] /= norm_e1;

    std::array<real_t, 3> e2 = {
        n_norm[1]*e1[2] - n_norm[2]*e1[1],
        n_norm[2]*e1[0] - n_norm[0]*e1[2],
        n_norm[0]*e1[1] - n_norm[1]*e1[0]
    };

    std::array<real_t, 3> d = {x - P0[0], y - P0[1], z - P0[2]};

    s = d[0]*n_norm[0] + d[1]*n_norm[1] + d[2]*n_norm[2]; //s

    std::array<real_t, 3> d_trans;
    for (int i = 0; i < 3; ++i)
        d_trans[i] = d[i] - s * n_norm[i];

    real_t r_cos_theta = d_trans[0]*e1[0] + d_trans[1]*e1[1] + d_trans[2]*e1[2];
    real_t r_sin_theta = d_trans[0]*e2[0] + d_trans[1]*e2[1] + d_trans[2]*e2[2];
    
    r = sqrt(r_cos_theta*r_cos_theta + r_sin_theta*r_sin_theta); //r
    theta = atan2(r_sin_theta, r_cos_theta); //theta
}




/**
 * Computes the deposition profile for a flux surface
 */
real_t NBIElectronHeatTerm::ComputeDepositionProfile(real_t ir, real_t ncold, real_t Tcold, real_t Ti, real_t ni, FVM::UnknownQuantityHandler* unknowns
){

    len_t ntheta = 100;
    len_t nphi = 100;
    real_t dtheta = 2.0*M_PI /ntheta;
    real_t dphi = 2.0*M_PI/nphi;
    
    real_t total_deposition = 0.0;
    real_t volume_element = 0.0;

    for (len_t itheta = 0; itheta < ntheta; itheta++){ //integrate by theta
        real_t theta = itheta*dtheta;
        for (len_t iphi = 0; iphi < nphi; iphi++){ //integrate by phi
        //For every point on that flux surface
            real_t phi = iphi*dphi;

            real_t R = radialGrid->GetFluxSurfaceR(ir, theta);
            real_t Z = radialGrid->GetFluxSurfaceZ(ir, theta);

            real_t x = R * cos(phi); 
            real_t y = R * sin(phi);
            real_t z = Z;

            // Compute the beamline coordinates r_B,theta_B,s_B for this/every point, so now we are in the beam system
            real_t r_B, theta_B, s_B;
            CartesianToCylindrical(x, y, z, P0, n, a, r_B, theta_B, s_B); 
            //Check if our flux surface point is within the beam
            if (s_B < 0 || s_B > s_max)
               continue;

            if (r_B > r_beam)
                continue;

            real_t lambda_s = ComputeMeanFreePath(ir,ncold,Tcold,Ti,ni); // At current flux surface point
            real_t survivalProb = ComputeSurvivalProbability(ir,s_B, unknowns); //Up to current flux surface point

            // Compute deposition probability density H(r, theta, phi)
            real_t H_r_theta_phi = (j_B / I_B) * (1.0 / lambda_s) * survivalProb;
        

            real_t J = ComputeConfigurationSpaceJacobian(ir, theta); 

            // Integrate over poloidal and toroidal angles
            total_deposition += H_r_theta_phi * J * dtheta * dphi; 
            volume_element += J * dtheta * dphi; 

        }
    }
    return total_deposition / volume_element;
}

/**
 * Set the Jacobian elements corresponding to this term.
 */
bool NBIElectronHeatTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId,
    FVM::Matrix *jac, const real_t *x,
    FVM::UnknownQuantityHandler* unknowns
) {
    for (len_t ir = 0; ir < nr; ir++) {
        auto [dP_dni, dP_dTi, dP_dTe, dP_dne] = Compute_dP_derivative(ir, unknowns);

        real_t dP_dX = 0.0;

        if (derivId == id_ion_density) //TODO: HOW DO I KNOW THESE ARE THE CORRECT IDS?
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
    for (len_t ir = 0; ir < nr; ir++)
        rhs[ir] += NBIHeatTerm[ir];
}

/**
 * Sets the matrix elements for this term.
 */

void NBIElectronHeatTerm::SetMatrixElements(FVM::Matrix *mat, real_t *rhs) {
    //TODO: I FELL LIKE THIS IS WRONG
     for (len_t ir = 0; ir < nr; ir++) {
        rhs[ir] += NBIHeatTerm[ir];  // Add heating to RHS
    }
}


}


















