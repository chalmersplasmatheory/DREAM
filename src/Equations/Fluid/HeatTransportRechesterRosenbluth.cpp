/**
 * Implementation of a Rechester-Rosenbluth operator for heat transport.
 */

#include "DREAM/Equations/Fluid/HeatTransportRechesterRosenbluth.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
HeatTransportRechesterRosenbluth::HeatTransportRechesterRosenbluth(
    FVM::Grid *grid, enum OptionConstants::momentumgrid_type mgtype,
    FVM::Interpolator1D *dB_B, FVM::UnknownQuantityHandler *unknowns
) : FVM::DiffusionTerm(grid), mgtype(mgtype), deltaBOverB(dB_B) {

    this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_T_cold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    
    AddUnknownForJacobian(unknowns, this->id_n_cold);
    AddUnknownForJacobian(unknowns, this->id_T_cold);
}

/**
 * Destructor.
 */
HeatTransportRechesterRosenbluth::~HeatTransportRechesterRosenbluth() {
    delete this->deltaBOverB;
}


/**
 * Rebuild the coefficients for this equation term.
 */
void HeatTransportRechesterRosenbluth::Rebuild(
    const real_t t, const real_t, FVM::UnknownQuantityHandler *unknowns
) {
    const real_t *dB_B = this->deltaBOverB->Eval(t);
    const len_t nr = this->grid->GetNr();
    const real_t mc2 = Constants::mc2inEV;

    const real_t *ncold = unknowns->GetUnknownData(this->id_n_cold);
    const real_t *Tcold = unknowns->GetUnknownData(this->id_T_cold);

    real_t R0 = this->grid->GetRadialGrid()->GetR0();
    if (isinf(R0))
        R0 = 1;

    const real_t PREFAC = 3.0/2.0 * sqrt(2.0*M_PI) * Constants::c * R0 / mc2;

    for (len_t ir = 0; ir < nr; ir++) {
        auto rg = this->grid->GetRadialGrid();

        const real_t B_Bmin = rg->GetFSA_B(ir);
        const real_t xiT0   = sqrt(1 - rg->GetBmin() / rg->GetBmax());
        const real_t q = 1.0;       // TODO (safety factor)
        const real_t Theta = Tcold[ir] / mc2;

        real_t D = PREFAC * q * dB_B[ir]*dB_B[ir] 
            * B_Bmin * (1-xiT0) * ncold[ir] * sqrt(Theta) * (1 - 5.0/8.0*Theta) / mc2;
        Drr(ir, 0, 0) += D;
    }
}

/**
 * Set jacobian of diffusion coefficients for this diffusion term.
 *
 * derivId:    ID of the quantity with respect to which the coefficient should
 *             be differentiated.
 * nMultiples: (not used).
 */
void HeatTransportRechesterRosenbluth::SetPartialDiffusionTerm(
    len_t derivId, len_t
) {
    ResetDifferentiationCoefficients();
}

