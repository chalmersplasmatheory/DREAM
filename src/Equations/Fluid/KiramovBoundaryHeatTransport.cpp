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
    : FVM::AdvectionTerm(g, true) {

    SetName("KiramovBoundaryHeatTransport");

    this->id_T_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    
    for(len_t i = 0; i < ions->GetNZ(); i++){
        this->Z = ions->GetZ(i);
        if (Z == 1) {
            this->mi = ions->GetIonSpeciesMass(i);
            break;
        }
    }
}


/**
 * Build the coefficients of this advection term.
 */
void KiramovBoundaryHeatTransport::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *x) {
    real_t *T_cold_term = x->GetUnknownData(id_T_cold);
    real_t *n_cold_term = x->GetUnknownData(id_n_cold);
          
    real_t T = (1-deltaRadialFlux[nr]) * T_cold_term[nr-1] * Constants::ec;
    real_t n = (1-deltaRadialFlux[nr]) * n_cold_term[nr-1];

    real_t c_s = sqrt(adb_index * Z * T / mi); // Ions sound speed 

	// NOTE: Here we DON'T want to add to the coefficient, as we would
	// normally want. There is AdvectionDiffusion term corresponding to
	// this B.C. which would be responsible for re-setting the coefficient,
	// and so we can manually reset it by doing this.
	// (adding is only needed when we have multiple terms writing to the
	// coefficient, but for a boundary condition that should not be the
	// case as we should only have 1 on BC in the radial direction)
    Fr(nr,0,0) = gamma * n * c_s * Constants::ec; 

    T_cold = T_cold_term;
    n_cold = n_cold_term;    
}

// Set derivatives of the advection term with respect to different unknowns quantities 
void KiramovBoundaryHeatTransport::SetPartialAdvectionTerm(len_t derivId, len_t /*nMultiples*/){
    ResetDifferentiationCoefficients();
    
    real_t T = (1-deltaRadialFlux[nr]) * T_cold[nr-1]* Constants::ec;
    real_t n = (1-deltaRadialFlux[nr]) * n_cold[nr-1];

    real_t c_s = sqrt(adb_index * Z * T / mi);

    if (derivId == id_T_cold) {
        dFr(nr,0,0,0) += n * gamma * sqrt(adb_index *Z /mi) * 0.5 * Constants::ec / sqrt(T) ;
    } else if (derivId == id_n_cold) {
        dFr(nr,0,0,0) += gamma * c_s * Constants::ec;
    }
}

