/**
 * Implementation of the Dreicer runaway rate.
 */

#include "DREAM/Equations/Fluid/DreicerRateTerm.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "DREAM/DREAMException.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
DreicerRateTerm::DreicerRateTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *uqn,
    RunawayFluid *rf, IonHandler *ions, enum dreicer_type type, real_t scaleFactor
) : EquationTerm(g), REFluid(rf), ions(ions), type(type),
    scaleFactor(scaleFactor) {

    this->AllocateGamma();

    this->id_E_field = uqn->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    this->id_n_cold  = uqn->GetUnknownID(OptionConstants::UQTY_N_COLD);
}

/**
 * Destructor.
 */
DreicerRateTerm::~DreicerRateTerm() {
    DeallocateGamma();
}

/**
 * Allocate memory for the runaway rate.
 */
void DreicerRateTerm::AllocateGamma() {
    this->gamma = new real_t[this->grid->GetNr()];
}

/**
 * Free memory for the runaway rate.
 */
void DreicerRateTerm::DeallocateGamma() {
    delete [] this->gamma;
}

/**
 * Method called when the grid has been rebuilt.
 */
bool DreicerRateTerm::GridRebuilt() {
    DeallocateGamma();
    AllocateGamma();

    return true;
}

/**
 * Evaluate the Connor-Hastie runaway rate.
 *
 * ir:   Index of radius for which to evaluate the runaway rate.
 * E:    Local electric field strength.
 * Zeff: Local effective plasma charge.
 */
real_t DreicerRateTerm::ConnorHastie(
    const len_t ir, const real_t E, const real_t ne, const real_t Zeff
) {
    if (E == 0)
        return 0;

    real_t C  = 1;

    real_t Ec = REFluid->GetConnorHastieField_COMPLETESCREENING(ir);
    real_t ED = REFluid->GetDreicerElectricField(ir);

    if (Ec >= E)
        return 0;

    real_t EEc   = E / Ec;
    real_t EED   = E / ED;

    real_t tauEE = REFluid->GetElectronCollisionTimeThermal(ir);

    real_t h    = 1.0/(3*(EEc-1)) * (EEc + 2*(EEc-2)*sqrt(EEc/(EEc-1)) - (Zeff-7)/(Zeff+1));
    real_t etaf = (M_PI/2 - asin(1-2/EEc));
    real_t eta  = EEc*EEc/(4*(EEc-1)) * etaf*etaf;
    real_t lmbd = 8*EEc*EEc*(1 - 1.0/(2*EEc) - sqrt(1-1/EEc));

    return C*ne/tauEE * pow(EED, 3/16*(1+Zeff)*h)
        * exp(-lmbd/(4*EED) - sqrt(eta*(1+Zeff)/EED));
}

/**
 * Calculate the Dreicer runaway rate at the current time.
 */
void DreicerRateTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *uqn) {
    real_t *E  = uqn->GetUnknownData(this->id_E_field);
    real_t *ne = uqn->GetUnknownData(this->id_n_cold);

    const len_t nr = this->grid->GetNr();

    for (len_t ir = 0; ir < nr; ir++) {
        real_t Zeff = ions->evaluateZeff(ir);

        switch (this->type) {
            case CONNOR_HASTIE:
                gamma[ir] = ConnorHastie(ir, E[ir], ne[ir], Zeff);
                break;

            case NEURAL_NETWORK:
                throw NotImplementedException("The Dreicer neural network has not been implemented yet.");

            default:
                throw DREAMException(
                    "Unrecognized Dreicer generation type: %d.",
                    this->type
                );
        }
    }
}

/**
 * Set the Jacobian elements corresponding to this term.
 */
void DreicerRateTerm::SetJacobianBlock(
    const len_t /*uqtyId*/, const len_t /*derivId*/, FVM::Matrix* /*jac*/, const real_t* /*qty*/
) {
    // TODO
}

/**
 * Set the linear operator matrix elements corresponding
 * to this term.
 */
void DreicerRateTerm::SetMatrixElements(FVM::Matrix*, real_t *rhs) {
    this->SetVectorElements(rhs, nullptr);
}

/**
 * Set the non-linear function vector for this term.
 */
void DreicerRateTerm::SetVectorElements(real_t *vec, const real_t*) {
    const len_t nr  = this->grid->GetNr();
    len_t offset = 0;
    
    for (len_t ir = 0; ir < nr; ir++) {
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();
        const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();

        // Insert at p=0, xi=1
        vec[offset + np1*(np2-1) + 0] = scaleFactor*this->gamma[ir];

        offset += np1*np2;
    }
}

