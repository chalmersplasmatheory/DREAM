#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/Fluid/KiramovBoundaryHeatTransport.hpp"

using namespace DREAM;

/**
 * Implementation of the Kiramov boundary condition for the heat transport
 * equation. This boundary condition is given by
 *
 *   q|| = -gamma * n_e * c_s(T_e) * T_e 
 *
 * where q|| is the parallel heat flux, c_s is the sound speed and gamma is the 
 * heat transmission coefficient
 */


/**
 * Constructor.
 */
KiramovBoundaryHeatTransport::KiramovBoundaryHeatTransport(FVM::Grid *g, FVM::UnknownQuantityHandler *unknowns, IonHandler *ions)
    : FVM::AdvectionTerm(g, true), unknowns(unknowns) {

    SetName("KiramovBoundaryHeatTransport");

    this->id_T_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    if(unknowns->HasUnknown(OptionConstants::UQTY_NI_DENS) && unknowns->HasUnknown(OptionConstants::UQTY_WI_ENER)){
        this->id_N_i = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
        this->id_W_i = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
        has_Ti = true;
    }
    this->id_jtot  = unknowns->GetUnknownID(OptionConstants::UQTY_J_TOT);



    // Loop through all ions to get the one with the minimum Z
    for (len_t i = 0; i < ions->GetNZ(); i++) {
        int currentZ = ions->GetZ(i); 
        if (currentZ < minZ) {        
            minZ = currentZ;          
            minIndex = i;             
        }
    }

    // Check if a valid minimum Z was found and set the ion quantities
    if (minIndex != -1) {
        this-> mi = ions->GetIonSpeciesMass(minIndex); 
        this-> Z = minZ;  

    }
    const real_t
        *Vp_fr = this->grid->GetVp_fr(this->grid->GetNr()),
        *Vp    = this->grid->GetVp(this->grid->GetNr()-1),
        dr    = this->grid->GetRadialGrid()->GetDr(this->grid->GetNr()-1);
    S_wo_coeff = Vp_fr[0] / (Vp[0] * dr);
}


/**
 * Build the coefficients of this advection term.
 */
void KiramovBoundaryHeatTransport::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *) {// TODO
    real_t *T_cold = unknowns->GetUnknownData(id_T_cold); 
    real_t *n_cold = unknowns->GetUnknownData(id_n_cold); 

    real_t T_i;
    if (has_Ti) {
        real_t *N_i = unknowns->GetUnknownData(id_N_i); 
        real_t *W_i = unknowns->GetUnknownData(id_W_i); 
        T_i = 2. / 3. * W_i[nr-1] / N_i[nr-1]; 
    } else 
        T_i = 0.;

    real_t T_e = T_cold[nr-1] * Constants::ec; // (1-deltaRadialFlux[nr]) * T_cold[nr-1] * Constants::ec;
    real_t n_e = n_cold[nr-1];                 // (1-deltaRadialFlux[nr]) * n_cold[nr-1];

    real_t c_s = sqrt((T_e + gamma * T_i) / mi); //sqrt(adb_index * Z * T / mi); // Ions sound speed 

    const real_t *jtot = this->unknowns->GetUnknownData(id_jtot);
    real_t mu0Ip = Constants::mu0 * TotalPlasmaCurrentFromJTot::EvaluateIpInsideR(this->grid->GetNr()-1,this->grid->GetRadialGrid(),jtot);
    real_t qR0 = this->grid->GetRadialGrid()->SafetyFactorNormalized(this->grid->GetNr()-1,mu0Ip);
    real_t L_par = 2 * M_PI * qR0;

	// NOTE: Here we DON'T want to add to the coefficient, as we would
	// normally want. There is AdvectionDiffusion term corresponding to
	// this B.C. which would be responsible for re-setting the coefficient,
	// and so we can manually reset it by doing this.
	// (adding is only needed when we have multiple terms writing to the
	// coefficient, but for a boundary condition that should not be the
	// case as we should only have 1 on BC in the radial direction)
    real_t q_par = 2 * kappa * n_e * c_s * Constants::ec;

    Fr(nr,0,0) = 2. / 3. * q_par / L_par * 1. / S_wo_coeff; 
}

// Set derivatives of the advection term with respect to different unknowns quantities 
void KiramovBoundaryHeatTransport::SetPartialAdvectionTerm(len_t derivId, len_t /*nMultiples*/){ 
    /* TODO: Deriv for j_tot? And is this ok for T_cold? We apply to T_cold...
    */
    ResetDifferentiationCoefficients();
    
    real_t *T_cold = unknowns->GetUnknownData(id_T_cold); 
    
    real_t T_e = T_cold[nr-1] * Constants::ec; // (1-deltaRadialFlux[nr]) * T_cold[nr-1]* Constants::ec;
    real_t n_e = n_cold[nr-1];                 // (1-deltaRadialFlux[nr]) * n_cold[nr-1];

    const real_t *jtot = this->unknowns->GetUnknownData(id_jtot);
    real_t mu0Ip = Constants::mu0 * TotalPlasmaCurrentFromJTot::EvaluateIpInsideR(this->grid->GetNr()-1,this->grid->GetRadialGrid(),jtot);
    real_t qR0 = this->grid->GetRadialGrid()->SafetyFactorNormalized(this->grid->GetNr()-1,mu0Ip);
    real_t L_par = 2 * M_PI * qR0;
    
    real_t T_i;
    real_t c_s;
    if (has_Ti) {
        real_t *N_i = unknowns->GetUnknownData(id_N_i); 
        real_t *W_i = unknowns->GetUnknownData(id_W_i); 
        T_i = 2. / 3. * W_i[nr-1] / N_i[nr-1]; 

        c_s = sqrt((T_e + gamma * T_i) / mi);
        if (derivId == id_W_i) {
            real_t d_c_s = 1./ 2. * 1. / sqrt((T_e + gamma * T_i) / mi) * gamma / mi * 2. / 3. / N_i[nr-1]; 
            real_t d_q_par = 2 * kappa * n_e * d_c_s * Constants::ec;
            dFr(nr,0,0,0) += 2. / 3. * d_q_par / L_par * 1. / S_wo_coeff; 
        } else if (derivId == id_N_i) {
            real_t d_c_s = 1./ 2. * 1. / sqrt((T_e + gamma * T_i) / mi) * gamma / mi * 2. / 3. * (- W_i[nr-1] / (N_i[nr-1] * N_i[nr-1])); 
            real_t d_q_par = 2 * kappa * n_e * d_c_s * Constants::ec;
            dFr(nr,0,0,0) += 2. / 3. * d_q_par / L_par * 1. / S_wo_coeff; 
        }
    } else 
        T_i = 0.;

    c_s = sqrt((T_e + gamma * T_i) / mi); //sqrt(adb_index * Z * T / mi); // Ions sound speed 

    if (derivId == id_T_cold) {
        real_t d_c_s = 1./ 2. * 1. / sqrt((T_e + gamma * T_i) / mi) * Constants::ec / mi;
        real_t d_q_par = 2 * kappa * n_e * d_c_s * Constants::ec;
        dFr(nr,0,0,0) += 2. / 3. * d_q_par / L_par * 1. / S_wo_coeff; //n * gamma * sqrt(adb_index *Z /mi) * 0.5 * Constants::ec / sqrt(T) ;
    } else if (derivId == id_n_cold) {
        real_t d_q_par = 2 * kappa * c_s * Constants::ec;
        dFr(nr,0,0,0) += 2. / 3. * d_q_par / L_par * 1. / S_wo_coeff; // gamma * c_s * Constants::ec;
    }
}
