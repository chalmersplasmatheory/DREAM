/**
 * Implementation of a Rechester-Rosenbluth operator for heat transport.
 */

#include <cmath>
#include "DREAM/Constants.hpp"
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

    SetName("HeatTransportRechesterRosenbluth");

    this->unknowns = unknowns;
    this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_T_cold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    
    AddUnknownForJacobian(unknowns, this->id_n_cold);
    AddUnknownForJacobian(unknowns, this->id_T_cold);

    AllocateDiffCoeff();
}

/**
 * Destructor.
 */
HeatTransportRechesterRosenbluth::~HeatTransportRechesterRosenbluth() {
    delete this->deltaBOverB;
    delete [] this->dD;
}


/**
 * Allocate memory for the differentiation coefficient.
 */
void HeatTransportRechesterRosenbluth::AllocateDiffCoeff() {
    const len_t nr = this->grid->GetNr();
    this->dD = new real_t[nr+1];
}

/**
 * Called whenever the grid is rebuilt.
 */
bool HeatTransportRechesterRosenbluth::GridRebuilt() {
    this->FVM::DiffusionTerm::GridRebuilt();

    delete [] this->dD;
    AllocateDiffCoeff();
    
    return true;
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

    FVM::RadialGrid *rg = this->grid->GetRadialGrid();

    const real_t PREFAC = 3.0 * sqrt(2.0*M_PI) * Constants::c * Constants::ec;
    for (len_t ir = 0; ir < nr+1; ir++) {
        real_t T=0, n=0;
        if(ir<nr){
            T += deltaRadialFlux[ir] * Tcold[ir];
            n += deltaRadialFlux[ir] * ncold[ir];
        } 
        if(ir>0){
            T += (1-deltaRadialFlux[ir]) * Tcold[ir-1];
            n += (1-deltaRadialFlux[ir]) * ncold[ir-1];
        }
        real_t Theta = T / mc2;
        
        const real_t B_Bmin = rg->GetFSA_B_f(ir);
        const real_t xiT0   = rg->GetXi0TrappedBoundary_fr(ir);
        real_t qR0;
        const real_t R0 = rg->GetR0();
        if(isinf(R0))
            qR0 = 1;       // TODO (safety factor)
        else
            qR0 = R0;       // TODO (safety factor)

        real_t D = PREFAC * dB_B[ir]*dB_B[ir] * B_Bmin * (1-xiT0*xiT0); // mc2;
        this->dD[ir] = D;
        
        Drr(ir, 0, 0) += qR0 * D * n * sqrt(Theta) * (1 - 5.0/8.0*Theta);
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

    const len_t nr = this->grid->GetNr();
    const real_t mc2 = Constants::mc2inEV;

    const real_t *ncold = unknowns->GetUnknownData(this->id_n_cold);
    const real_t *Tcold = unknowns->GetUnknownData(this->id_T_cold);

    for (len_t ir = 0; ir < nr+1; ir++) {
        real_t T=0, n=0;
        if(ir<nr){
            T += deltaRadialFlux[ir] * Tcold[ir];
            n += deltaRadialFlux[ir] * ncold[ir];
        } 
        if(ir>0){
            T += (1-deltaRadialFlux[ir]) * Tcold[ir-1];
            n += (1-deltaRadialFlux[ir]) * ncold[ir-1];
        }
        real_t Theta = T / mc2;
        real_t qR0;
        const real_t R0 = this->grid->GetRadialGrid()->GetR0();
        if(isinf(R0))
            qR0 = 1;       // TODO: safety factor with j_tot jacobian
        else
            qR0 = R0;       // TODO: safety factor with j_tot jacobian

        if (derivId == this->id_n_cold)
            dDrr(ir, 0, 0) = qR0 * this->dD[ir] * sqrt(Theta) * (1 - 5.0/8.0*Theta);
        else if (derivId == this->id_T_cold)
            dDrr(ir, 0, 0) = qR0 * this->dD[ir]/mc2 * n * (0.5 - 15.0/16.0*Theta) / sqrt(Theta);
    }
}
