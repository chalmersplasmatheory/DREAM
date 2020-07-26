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
    this->id_T_cold  = uqn->GetUnknownID(OptionConstants::UQTY_T_COLD);
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
    this->gamma           = new real_t[this->grid->GetNr()];
    this->EED_dgamma_dEED = new real_t[this->grid->GetNr()];
}

/**
 * Free memory for the runaway rate.
 */
void DreicerRateTerm::DeallocateGamma() {
    delete [] this->EED_dgamma_dEED;
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
 * Calculate the Dreicer runaway rate at the current time.
 */
void DreicerRateTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *uqn) {
    const len_t nr = this->grid->GetNr();

    for (len_t ir = 0; ir < nr; ir++)
        this->gamma[ir] = REFluid->GetDreicerRunawayRate(ir);

    if (this->type == CONNOR_HASTIE) {
        ConnorHastie *ch = REFluid->GetConnorHastieRunawayRate();

        const real_t *E  = uqn->GetUnknownData(id_E_field);
        const real_t *n  = uqn->GetUnknownData(id_n_cold);

        for (len_t ir = 0; ir < nr; ir++) {
            real_t EED  = E[ir] / REFluid->GetDreicerElectricField(ir);
            real_t Zeff = this->ions->evaluateZeff(ir);

            this->EED_dgamma_dEED[ir] = EED * ch->Diff_EED(ir, E[ir], n[ir], Zeff);
        }
    } else if (this->type == NEURAL_NETWORK) {
        // TODO
    }
}

/**
 * Set the Jacobian elements corresponding to this term.
 */
void DreicerRateTerm::SetJacobianBlock(
    const len_t, const len_t derivId, FVM::Matrix *jac, const real_t *qty
) {
    const len_t nr = this->grid->GetNr();

    if (derivId == id_E_field || derivId == id_n_cold || derivId == id_T_cold) {
        real_t s = (derivId == id_n_cold ? -1.0 : 1.0);

        for (len_t ir = 0; ir < nr; ir++) {
            if (qty[ir] == 0) continue;

            real_t v = s*this->EED_dgamma_dEED[ir] / qty[ir];
            jac->SetElement(ir, ir, v);
        }
    }
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

