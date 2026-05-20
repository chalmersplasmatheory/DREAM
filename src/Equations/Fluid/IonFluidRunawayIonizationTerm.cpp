/**
 * Implementation of the kinetic ionization term, assuming
 * a distribution function ~delta(p - p*), for some given p*,
 *
 * Note that this equation is applied to a single _ion species_,
 * (and to all its charge states).
 */

#include "DREAM/Equations/Fluid/IonFluidRunawayIonizationTerm.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
IonFluidRunawayIonizationTerm::IonFluidRunawayIonizationTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ihdl, const len_t iIon, real_t scaleFactor=1.0
) : IonEquationTerm<FVM::EquationTerm>(g, ihdl, iIon), u(u), scaleFactor(scaleFactor) {

    SetName("IonFluidRunawayIonizationTerm");

    id_ions = u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    id_nre = u->GetUnknownID(OptionConstants::UQTY_N_RE);
    AddUnknownForJacobian(u, id_ions);
    AddUnknownForJacobian(u, id_nre);

    real_t p0 = CHARACTERISTIC_RUNAWAY_MOMENTUM;
    preFactor = Constants::c * p0 / sqrt(1 + p0*p0);

    nr = g->GetNr();
    weights = new real_t[(Zion+1)*nr];

    tableIndexIon = IonKineticIonizationTerm::GetTableIndex(Zion);
}

/**
 * Destructor.
 */
IonFluidRunawayIonizationTerm::~IonFluidRunawayIonizationTerm() {
    delete [] weights;
}

void IonFluidRunawayIonizationTerm::Rebuild(
    const real_t /* t */, const real_t /* dt */, FVM::UnknownQuantityHandler* /* u */
) {
    // set weights
    for (len_t iZ = 0; iZ < Zion; iZ++) {
        const real_t *params = IonKineticIonizationTerm::kinetic_rate_table[tableIndexIon].params + iZ *  IonKineticIonizationTerm::nParamsForFit;
        for (len_t ir = 0; ir < nr; ir++) {
            weights[ir*(Zion+1)+iZ] = scaleFactor * preFactor * IonKineticIonizationTerm::EvaluateIonizationCrossSection(CHARACTERISTIC_RUNAWAY_MOMENTUM, params);
        }
    }
    for (len_t ir = 0; ir < nr; ir++)
        weights[ir*(Zion+1)+Zion] = 0;  // fully ionized ions do not ionize further...
}


bool IonFluidRunawayIonizationTerm::SetCSJacobianBlock(
    const len_t /* uqtyId */, const len_t derivId, FVM::Matrix *jac, const real_t* /* nions */,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
    if ( !HasJacobianContribution(derivId) )
        return false;

    if (derivId == id_ions)
        this->SetCSMatrixElements(jac, nullptr, iIon, Z0, rOffset);
    if (derivId == id_nre)
        for (len_t ir = 0; ir < nr; ir++){
            real_t nion = ions->GetIonDensity(ir, iIon, Z0);
            jac->SetElement(rOffset+ir, ir, -weights[ir*(Zion+1)+Z0] * nion);
            if (Z0 > 0) {
                nion = ions->GetIonDensity(ir, iIon, Z0-1);
                jac->SetElement(rOffset+ir, ir, weights[ir*(Zion+1)+Z0-1] * nion);
            }
        }
    return true;
}


/**
 * Sets the matrix elements of this equation term
 */
void IonFluidRunawayIonizationTerm::SetCSMatrixElements(
    FVM::Matrix *mat, real_t* /* rhs */, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
    real_t *nre = u->GetUnknownData(id_nre);
    for (len_t ir = 0; ir < nr; ir++){
        mat->SetElement(rOffset+ir, rOffset+ir, -weights[ir*(Zion+1)+Z0] * nre[ir]);
        if (Z0 > 0)
            mat->SetElement(rOffset+ir, rOffset+ir-nr, weights[ir*(Zion+1)+Z0-1] * nre[ir]);
    }
}


/**
 * Sets vector elements for this ion and charge state
 */
void IonFluidRunawayIonizationTerm::SetCSVectorElements(
    real_t *vec, const real_t *nions, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
    real_t *nre = u->GetUnknownData(id_nre);
    for(len_t ir = 0; ir < nr; ir++){
        vec[rOffset+ir] -= weights[ir*(Zion+1)+Z0] * nre[ir] * nions[rOffset+ir];
        if (Z0 > 0)
            vec[rOffset+ir] += weights[ir*(Zion+1)+Z0-1] * nre[ir] * nions[rOffset+ir-nr];
    }
}
