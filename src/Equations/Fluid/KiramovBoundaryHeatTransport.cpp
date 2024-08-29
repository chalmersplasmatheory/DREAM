#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/Fluid/KiramovBoundaryHeatTransport.hpp"

using namespace DREAM;

/**
 * Implementation of the Kiramov boundary condition for the heat transport
 * equation. This boundary condition is given by
 *
 *   q|| = -g n_e c_s(T_e) * T_e 
 *
 * where q|| is the parallel heat flux, c_s is the sound speed and gamma is the 
 * heat transmission coefficient
 */


/**
 * Constructor.
 */
KiramovBoundaryHeatTransport::KiramovBoundaryHeatTransport(FVM::Grid *g, FVM::UnknownQuantityHandler *unknowns)
    : FVM::AdvectionTerm(g, true) {

    SetName("KiramovBoundaryHeatTransport");

    this->id_T_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);

}


/**
 * Build the coefficients of this advection term.
 */
void KiramovBoundaryHeatTransport::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *x) {
    real_t *T_cold_term = x->GetUnknownData(id_T_cold);
    real_t *n_cold_term = x->GetUnknownData(id_n_cold);
          
    real_t T = (1-deltaRadialFlux[nr]) * T_cold_term[nr-1];
    real_t n = (1-deltaRadialFlux[nr]) * n_cold_term[nr-1];

    real_t c_s = sqrt(adb_index * T / Constants::me);

    Fr(nr,0,0) += gamma * n * c_s * T; 

    T_cold = T_cold_term;
    n_cold = n_cold_term;    
}

// Set derivatives of the advection term with respect to different unknowns quantities 
void KiramovBoundaryHeatTransport::SetPartialAdvectionTerm(len_t derivId, len_t /*nMultiples*/){
    ResetDifferentiationCoefficients();
    
    real_t T = (1-deltaRadialFlux[nr]) * T_cold[nr-1];
    real_t n = (1-deltaRadialFlux[nr]) * n_cold[nr-1];

    real_t c_s = sqrt(adb_index * T / Constants::me);

    if(derivId == id_T_cold){
        dFr(nr,0,0,0) += n * gamma * sqrt(adb_index/Constants::me) * 0.5 * T / sqrt(T) ;
        
    } else if(derivId == id_n_cold){
        dFr(nr,0,0,0) += gamma * c_s * T;
    }
}

