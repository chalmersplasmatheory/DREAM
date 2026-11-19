/**
 * Implementation of the thermal energy that atoms carry between the neutral
 * heat W_n and the ion heat W_i when they are ionized (Z0=0 -> 1) or
 * recombine (Z0=1 -> 0):
 */
#include "DREAM/Equations/Fluid/NeutralIonizationEnergyTransferTerm.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonKineticIonizationTerm.hpp"
#include "DREAM/Equations/Fluid/IonKineticIonizationTerm.hpp"
#include <cmath>
#include <limits>

using namespace DREAM;

/**
 * Constructor.
 */
NeutralIonizationEnergyTransferTerm::NeutralIonizationEnergyTransferTerm(
    FVM::Grid *grid, len_t iz, FVM::UnknownQuantityHandler *unknowns,
    IonHandler *ions, RateHandler *ratehandler, bool neutralEquation,
    bool includeFluidREIonization
) : FVM::EquationTerm(grid), unknowns(unknowns), ions(ions),
    ratehandler(ratehandler), iz(iz), neutralEquation(neutralEquation) {

    SetName("NeutralIonizationEnergyTransferTerm");

    this->id_ions = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this->id_ne   = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Te   = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_n_re = unknowns->GetUnknownID(OptionConstants::UQTY_N_RE);
    this->id_Wn   = unknowns->GetUnknownID(OptionConstants::UQTY_WN_ENER);
    this->id_Nn   = unknowns->GetUnknownID(OptionConstants::UQTY_NN_DENS);
    this->id_Wi   = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_Ni   = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);

    AddUnknownForJacobian(unknowns, id_ions);
    AddUnknownForJacobian(unknowns, id_ne);
    AddUnknownForJacobian(unknowns, id_Te);
    AddUnknownForJacobian(unknowns, id_n_re);
    AddUnknownForJacobian(unknowns, id_Wn);
    AddUnknownForJacobian(unknowns, id_Nn);
    AddUnknownForJacobian(unknowns, id_Wi);
    AddUnknownForJacobian(unknowns, id_Ni);

    this->reIonizationRate = 0;

    if (includeFluidREIonization){
        const real_t p0 = 20.0;
        const len_t tableIndex = IonKineticIonizationTerm::GetTableIndex(ions->GetZ(iz));
         const real_t *params = IonKineticIonizationTerm::kinetic_rate_table[tableIndex].params + Z0 * IonKineticIonizationTerm::nParamsForFit;
        reIonizationRate = reionizScale * Constants::c * p0 / sqrt(1 + p0*p0)
        * IonKineticIonizationTerm::EvaluateIonizationCrossSection(p0, params);
    }

    const len_t nr = grid->GetNr();
    this->Ion   = new real_t[nr];
    this->Ion_N = new real_t[nr];
    this->Ion_T = new real_t[nr];
    this->Rec   = new real_t[nr];
    this->Rec_N = new real_t[nr];
    this->Rec_T = new real_t[nr];
}

/**
 * Destructor.
 */
NeutralIonizationEnergyTransferTerm::~NeutralIonizationEnergyTransferTerm() {
    delete [] this->Ion;
    delete [] this->Ion_N;
    delete [] this->Ion_T;
    delete [] this->Rec;
    delete [] this->Rec_N;
    delete [] this->Rec_T;
}

/**
 * Rebuild the rate coefficients.
 */
void NeutralIonizationEnergyTransferTerm::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler*
) {
    const len_t nr = grid->GetNr();
    const real_t *ne = unknowns->GetUnknownData(id_ne);
    const real_t *Te = unknowns->GetUnknownData(id_Te);

    for (len_t ir = 0; ir < nr; ir++) {
        Ion[ir]   = ratehandler->GetSCD(iz)->Eval(0,ne[ir], Te[ir]);
        Ion_N[ir] = (ratehandler->GetSCD(iz)->Eval_deriv_n(0,ne[ir], Te[ir]));
        Ion_T[ir] = (ratehandler->GetSCD(iz)->Eval_deriv_T(0,ne[ir], Te[ir]));

        Rec[ir]   = ratehandler->GetACD(iz)->Eval(1, ne[ir], Te[ir]);
        Rec_N[ir] = (ratehandler->GetACD(iz)->Eval_deriv_n(1, ne[ir], Te[ir]));
        Rec_T[ir] = (ratehandler->GetACD(iz)->Eval_deriv_T(1, ne[ir], Te[ir]));
    }
}

/**
 * Set the elements of the vector corresponding to this term.
 */
void NeutralIonizationEnergyTransferTerm::SetVectorElements(
    real_t *vec, const real_t*
) {
    const len_t nr = grid->GetNr();
    const real_t *ni = unknowns->GetUnknownData(id_ions);
    const real_t *ne = unknowns->GetUnknownData(id_ne);
    const real_t *nre = unknowns->GetUnknownData(id_n_re);
    const len_t offset0 = ions->GetIndex(iz, 0) * nr;
    const len_t offset1 = ions->GetIndex(iz, 1) * nr;
    const real_t c = (neutralEquation ? -1 : 1) * 1.5 * Constants::ec;

    for (len_t ir = 0; ir < nr; ir++) {
        const real_t Tn = ratehandler->ResolveSpeciesTemperature(iz, 0, ir);
        const real_t Ti = ratehandler->ResolveSpeciesTemperature(iz, 1, ir);

        const real_t ionization    = (Ion[ir] * ne[ir] + reIonizationRate * nre[ir]) * ni[offset0 + ir];
        const real_t recombination = Rec[ir] * ne[ir] * ni[offset1 + ir];

        vec[iz*nr + ir] += c * (ionization * Tn - recombination * Ti);
    }
}

/**
 * Build block of Jacobian matrix.
 *
 * derivId: ID of unknown quantity with respect to which differentiation
 *          should be carried out.
 * jac:     Jacobian matrix to build.
 */
bool NeutralIonizationEnergyTransferTerm::SetJacobianBlock(
    const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*
) {
    if (!HasJacobianContribution(derivId))
        return false;

    const len_t nr = grid->GetNr();
    const real_t *ni = unknowns->GetUnknownData(id_ions);
    const real_t *ne = unknowns->GetUnknownData(id_ne);
    const real_t *nre = unknowns->GetUnknownData(id_n_re);
    const len_t offset0 = ions->GetIndex(iz, 0) * nr;
    const len_t offset1 = ions->GetIndex(iz, 1) * nr;
    const real_t c = (neutralEquation ? -1 : 1) * 1.5 * Constants::ec;

    for (len_t ir = 0; ir < nr; ir++) {
        const len_t row = iz*nr + ir;
        const real_t n0 = ni[offset0 + ir];
        const real_t n1 = ni[offset1 + ir];
        const real_t Tn = ratehandler->ResolveSpeciesTemperature(iz, 0, ir);
        const real_t Ti = ratehandler->ResolveSpeciesTemperature(iz, 1, ir);

        // Ionization and recombination rates per atom [1/s]
        const real_t ionizationPerN0    = Ion[ir] * ne[ir] + reIonizationRate * nre[ir];
        const real_t recombinationPerN1 = Rec[ir] * ne[ir];

        if (derivId == id_ions) {
            jac->SetElement(row, offset0 + ir, c * ionizationPerN0 * Tn);
            jac->SetElement(row, offset1 + ir, -c * recombinationPerN1 * Ti);
        }
        else if (derivId == id_ne)
            jac->SetElement(row, ir, c * (
                (Ion[ir] + Ion_N[ir] * ne[ir]) * n0 * Tn
                - (Rec[ir] + Rec_N[ir] * ne[ir]) * n1 * Ti
            ));
        else if (derivId == id_Te)
            jac->SetElement(row, ir, c * ne[ir] * (Ion_T[ir] * n0 * Tn - Rec_T[ir] * n1 * Ti));
        else if (derivId == id_n_re)
            jac->SetElement(row, ir, c * reIonizationRate * n0 * Tn);

        //Adds the partial derivatives with respect to the species temperatures (Tn and Ti)
        ratehandler->SetSpeciesTemperatureJacobian(derivId, jac, row, ir, iz, 0, c * ionizationPerN0 * n0);
        ratehandler->SetSpeciesTemperatureJacobian(derivId, jac, row, ir, iz, 1, -c * recombinationPerN1 * n1);
    }

    return true;
}
