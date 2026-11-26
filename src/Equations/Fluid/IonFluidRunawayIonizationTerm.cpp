/**
 * Implementation of the kinetic ionization term, assuming
 * a distribution function ~delta(p - p*), for some given p*,
 *
 * Note that this equation is applied to a single _ion species_,
 * (and to all its charge states).
 */

#include "DREAM/Equations/Fluid/IonFluidRunawayIonizationTerm.hpp"
#include <cmath>

using namespace DREAM;

/**
 * Constructor.
 */
IonFluidRunawayIonizationTerm::IonFluidRunawayIonizationTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ihdl, const len_t iIon,
    real_t scaleFactor, const FVM::Grid *exteriorGrid
) : IonEquationTerm<FVM::EquationTerm>(g, ihdl, iIon), u(u), exteriorGrid(exteriorGrid), scaleFactor(scaleFactor) {

    SetName("IonFluidRunawayIonizationTerm");

    id_ions = u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    AddUnknownForJacobian(u, id_ions);
    if (exteriorGrid == nullptr)
        idRunaways = {u->GetUnknownID(OptionConstants::UQTY_N_RE)};
    else
        idRunaways = {
            u->GetUnknownID(OptionConstants::UQTY_N_EXTERNAL_POS),
            u->GetUnknownID(OptionConstants::UQTY_N_EXTERNAL_NEG)
        };
    for (len_t id : idRunaways)
        AddUnknownForJacobian(u, id);

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
    const real_t *params = IonKineticIonizationTerm::kinetic_rate_table[tableIndexIon].params;
    for (len_t ir = 0; ir < nr; ir++) {
        real_t p = CHARACTERISTIC_RUNAWAY_MOMENTUM;
        if (exteriorGrid != nullptr) {
            const auto *mg = exteriorGrid->GetMomentumGrid(ir);
            p = mg->GetP1_f(mg->GetNp1());
        }
        const real_t v = Constants::c * p / std::hypot(1.0,p);
        for (len_t Z0 = 0; Z0 < Zion; Z0++) {
            const real_t *fit = params + Z0*IonKineticIonizationTerm::nParamsForFit;
            weights[ir*(Zion+1)+Z0] = scaleFactor * v * IonKineticIonizationTerm::EvaluateIonizationCrossSection(p,fit);
        }
        weights[ir*(Zion+1)+Zion] = 0;
    }
}

real_t IonFluidRunawayIonizationTerm::GetRunawayDensity(len_t ir) const {
    real_t density = 0;
    for (len_t id : idRunaways)
        density += u->GetUnknownData(id)[ir];
    return density;
}

/** Write the per-target rate (s^-1) for this species into a charge-state array. */
void IonFluidRunawayIonizationTerm::SetRateCoefficients(real_t *coefficients) const {
    const len_t offset = ions->GetIndex(iIon,0)*nr;
    for (len_t ir = 0; ir < nr; ir++) {
        const real_t density = GetRunawayDensity(ir);
        for (len_t Z0 = 0; Z0 <= Zion; Z0++)
            coefficients[offset+Z0*nr+ir] = weights[ir*(Zion+1)+Z0] * density;
    }
}

bool IonFluidRunawayIonizationTerm::SetCSJacobianBlock(
    const len_t /* uqtyId */, const len_t derivId, FVM::Matrix *jac, const real_t* /* nions */,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
    if ( !HasJacobianContribution(derivId) )
        return false;

    if (derivId == id_ions)
        this->SetCSMatrixElements(jac, nullptr, iIon, Z0, rOffset);
    else
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
    for (len_t ir = 0; ir < nr; ir++){
        const real_t density = GetRunawayDensity(ir);
        mat->SetElement(rOffset+ir, rOffset+ir, -weights[ir*(Zion+1)+Z0] * density);
        if (Z0 > 0)
            mat->SetElement(rOffset+ir, rOffset+ir-nr, weights[ir*(Zion+1)+Z0-1] * density);
    }
}


/**
 * Sets vector elements for this ion and charge state
 */
void IonFluidRunawayIonizationTerm::SetCSVectorElements(
    real_t *vec, const real_t *nions, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
    for(len_t ir = 0; ir < nr; ir++){
        const real_t density = GetRunawayDensity(ir);
        vec[rOffset+ir] -= weights[ir*(Zion+1)+Z0] * density * nions[rOffset+ir];
        if (Z0 > 0)
            vec[rOffset+ir] += weights[ir*(Zion+1)+Z0-1] * density * nions[rOffset+ir-nr];
    }
}
