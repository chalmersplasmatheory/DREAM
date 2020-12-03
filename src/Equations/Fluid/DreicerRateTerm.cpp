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
) : EquationTerm(g), RunawaySourceTerm(g, uqn), unknowns(uqn), REFluid(rf), ions(ions), type(type),
    scaleFactor(scaleFactor) {

    this->AllocateGamma();

    this->id_E_field = uqn->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    this->id_n_cold  = uqn->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_n_tot   = uqn->GetUnknownID(OptionConstants::UQTY_N_TOT);
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

    if (this->type == CONNOR_HASTIE || this->type == CONNOR_HASTIE_NOCORR) {
        ConnorHastie *ch = REFluid->GetConnorHastieRunawayRate();

        const real_t *E  = uqn->GetUnknownData(id_E_field);
        const real_t *n  = uqn->GetUnknownData(id_n_cold);

        for (len_t ir = 0; ir < nr; ir++) {
            real_t EED  = E[ir] / REFluid->GetDreicerElectricField(ir);
            real_t Zeff = this->ions->GetZeff(ir);

            this->EED_dgamma_dEED[ir] = EED * ch->Diff_EED(ir, E[ir], n[ir], Zeff);
        }
    }/* else if (this->type == NEURAL_NETWORK) {
        DreicerNeuralNetwork *dnn = REFluid->GetDreicerNeuralNetwork();

        const real_t *E      = uqn->GetUnknownData(id_E_field);
        const real_t *ntot   = uqn->GetUnknownData(id_n_tot);
        const real_t *T_cold = uqn->GetUnknownData(id_T_cold);

        for (len_t ir = 0; ir < nr; ir++) {
            real_t EED  = E[ir] / REFluid->GetDreicerElectricField(ir);
            this->EED_dgamma_dEED[ir] = EED * dnn->Diff_EED(ir, E[ir], ntot[ir], T_cold[ir]);
        }
    }*/

    this->data_E_field = uqn->GetUnknownData(id_E_field);
    this->data_n_cold  = uqn->GetUnknownData(id_n_cold);
    this->data_n_tot   = uqn->GetUnknownData(id_n_tot);
    this->data_T_cold  = uqn->GetUnknownData(id_T_cold);
}

/**
 * Set the Jacobian elements corresponding to this term.
 */
void DreicerRateTerm::SetJacobianBlock(
    const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*
) {
    const len_t nr = this->grid->GetNr();

    if (type == NEURAL_NETWORK) {
        // Numerical derivative
        if (derivId == id_E_field || derivId == id_n_tot || derivId == id_T_cold) {
            const real_t h = 1e-3;

            for (len_t ir = 0; ir < nr; ir++) {
                DreicerNeuralNetwork *dnn = this->REFluid->GetDreicerNeuralNetwork();

                real_t g0 = this->gamma[ir], g1, v;
                if (derivId == id_E_field) {
                    v = data_E_field[ir]*h;
                    g1 = dnn->RunawayRate(ir, data_E_field[ir]+v, data_n_tot[ir], data_T_cold[ir]);
                } else if (derivId == id_n_tot) {
                    v  = data_n_tot[ir]*h;
                    g1 = dnn->RunawayRate(ir, data_E_field[ir], data_n_tot[ir]+v, data_T_cold[ir]);
                } else if (derivId == id_T_cold) {
                    v  = data_T_cold[ir]*h;
                    g1 = dnn->RunawayRate(ir, data_E_field[ir], data_n_tot[ir], data_T_cold[ir]+v);
                } else  /* No effect (due to earlier if) -- but keeps LLVM from complaining */
                    break;

                const len_t xiIndex = this->GetXiIndexForEDirection(ir);
                const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();
                real_t V = GetVolumeScaleFactor(ir);

                // Check if the quantity w.r.t. which we differentiate is a
                // fluid quantity, in which case it has np1=0, xiIndex=0
                len_t np1_op = np1, xiIndex_op = xiIndex;
                if (this->unknowns->GetUnknown(derivId)->NumberOfElements() == nr) {
                    np1_op = 1;
                    xiIndex_op = 0;
                }

                real_t dg;
                if (v == 0) dg = 0;
                else dg = (g1-g0)/v;

                // Place particles in p=0, xi=1
                jac->SetElement(ir + np1*xiIndex, ir + np1_op*xiIndex_op, this->scaleFactor * dg * V);
            }
        }
    } else {
        if (derivId == id_E_field || derivId == id_T_cold) {
            const real_t *data;

            if      (derivId == id_E_field) data = this->data_E_field;
            else if (derivId == id_T_cold)  data = this->data_T_cold;

            for (len_t ir = 0; ir < nr; ir++) {
                if (data[ir] == 0) continue;

                const len_t xiIndex = this->GetXiIndexForEDirection(ir);
                const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();
                real_t V = GetVolumeScaleFactor(ir);

                // Check if the quantity w.r.t. which we differentiate is a
                // fluid quantity, in which case it has np1=0, xiIndex=0
                len_t np1_op = np1, xiIndex_op = xiIndex;
                if (this->unknowns->GetUnknown(derivId)->NumberOfElements() == nr) {
                    np1_op = 1;
                    xiIndex_op = 0;
                }

                real_t v = this->EED_dgamma_dEED[ir] / data[ir];

                // Place particles in p=0, xi=1
                jac->SetElement(ir + np1*xiIndex, ir + np1_op*xiIndex_op, this->scaleFactor*v*V);
            }
        } else if (derivId == id_n_cold) {
            const real_t *n = this->data_n_cold;

            for (len_t ir = 0; ir < nr; ir++) {
                if (n[ir] == 0) continue;

                const len_t xiIndex = this->GetXiIndexForEDirection(ir);
                const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();
                real_t V = GetVolumeScaleFactor(ir);

                // Check if the quantity w.r.t. which we differentiate is a
                // fluid quantity, in which case it has np1=0, xiIndex=0
                len_t np1_op = np1, xiIndex_op = xiIndex;
                if (this->unknowns->GetUnknown(derivId)->NumberOfElements() == nr) {
                    np1_op = 1;
                    xiIndex_op = 0;
                }

                real_t v = (this->gamma[ir] - this->EED_dgamma_dEED[ir]) / n[ir];

                // Place particles in p=0, xi=1
                jac->SetElement(ir + np1*xiIndex, ir + np1_op*xiIndex_op, this->scaleFactor*v*V);
            }
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
        const len_t xiIndex = this->GetXiIndexForEDirection(ir);
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();
        const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();
        real_t V = GetVolumeScaleFactor(ir);

        // Insert at p=0, xi=+1 (or -1 if E is negative)
        vec[offset + np1*xiIndex + 0] += scaleFactor*this->gamma[ir] * V;

        offset += np1*np2;
    }
}

