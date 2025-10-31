/**
 * Implementation of a heat diffusion operator taking the form
 *
 *   d/dr ( V'*D*n * dT/dr )
 *
 * This operator should be applied to 'T_cold' (or 'T_hot').
 */

#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/HeatTransportDiffusion.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
HeatTransportDiffusion::HeatTransportDiffusion(
    FVM::Grid *grid, enum OptionConstants::momentumgrid_type mgtype,
    FVM::Interpolator1D *D, FVM::UnknownQuantityHandler *unknowns,
	const len_t id_n
) : FVM::DiffusionTerm(grid), mgtype(mgtype), coeffD(D) {

    SetName("HeatTransportDiffusion");

    this->unknowns = unknowns;
    this->id_n = id_n;
    
    AddUnknownForJacobian(unknowns, this->id_n);

    AllocateDiffCoeff();
}

/**
 * Destructor.
 */
HeatTransportDiffusion::~HeatTransportDiffusion() {
    delete this->coeffD;
    delete [] this->dD;
}


/**
 * Allocate memory for the differentiation coefficient.
 */
void HeatTransportDiffusion::AllocateDiffCoeff() {
    const len_t nr = this->grid->GetNr();
    this->dD = new real_t[nr+1];
}

/**
 * Called whenever the grid is rebuilt.
 */
bool HeatTransportDiffusion::GridRebuilt() {
    this->FVM::DiffusionTerm::GridRebuilt();

    delete [] this->dD;
    AllocateDiffCoeff();
    
    return true;
}

/**
 * Rebuild the coefficients for this equation term.
 */
void HeatTransportDiffusion::Rebuild(
    const real_t t, const real_t, FVM::UnknownQuantityHandler *unknowns
) {
    const real_t *D = this->coeffD->Eval(t);
    const len_t nr = this->grid->GetNr();

    const real_t *n = unknowns->GetUnknownData(this->id_n);

    for (len_t ir = 0; ir < nr+1; ir++) {
        real_t na=0;
        if(ir<nr)
            na += deltaRadialFlux[ir] * n[ir];
        if(ir>0)
            na += (1-deltaRadialFlux[ir]) * n[ir-1];

        // Factor ec (=elementary charge) to convert from
        // eV to joule
        this->dD[ir] = 1.5 * Constants::ec * D[ir];
        
        Drr(ir, 0, 0) += 1.5 * Constants::ec * D[ir] * na;
    }
}

/**
 * Set jacobian of diffusion coefficients for this diffusion term.
 *
 * derivId:    ID of the quantity with respect to which the coefficient should
 *             be differentiated.
 * nMultiples: (not used).
 */
void HeatTransportDiffusion::SetPartialDiffusionTerm(
    len_t derivId, len_t
) {
    if (derivId != this->id_n)
        return;

    ResetDifferentiationCoefficients();

    const len_t nr = this->grid->GetNr();
    for (len_t ir = 0; ir < nr+1; ir++)
        dDrr(ir, 0, 0) = this->dD[ir];
}

