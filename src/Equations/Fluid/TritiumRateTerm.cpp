#include "DREAM/Equations/Fluid/TritiumRateTerm.hpp"

/**
 * Implementation of a class which represents the Tritium runaway rate contribution to the n_re equation.
 * Employs the analytical generation rate calculated by RunawayFluid.
 */
using namespace DREAM;

/**
 * Constructor.
 *
 * g:           Grid on which this term lives.
 * uqn:         UnknownQuantityHandler object.
 * iIon:        Index of ion to which this term should be applied.
 * rf:          RunawayFluid object.
 * ions:        IonHandler object.
 * scaleFactor: Factor by which to rescale the equation term.
 */
TritiumRateTerm::TritiumRateTerm(
    FVM::Grid *g, IonHandler *ions, FVM::UnknownQuantityHandler *uqn, const len_t iIon,
    RunawayFluid *rf, real_t scaleFactor
) : IonEquationTerm<FVM::EquationTerm>(g, ions, iIon), RunawaySourceTerm(g,uqn),
    REFluid(rf), scaleFactor(scaleFactor) { }


/**
 * Set elements of the Jacobian corresponding to this equation term.
 */
void TritiumRateTerm::SetCSJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac,
    const real_t*,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
    if (uqtyId == derivId)
        this->SetCSMatrixElements(jac, nullptr, iIon, Z0, rOffset);
}

/**
 * Set elements of the linear operator matrix.
 */
void TritiumRateTerm::SetCSMatrixElements(
    FVM::Matrix *mat, real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
) {
    this->SetCSElements_internal(
        [&mat](const len_t i, const len_t j, const real_t v) { mat->SetElement(i, j, v); },
        iIon, Z0, rOffset
    );
}

void TritiumRateTerm::SetCSVectorElements(
    real_t *vec, const real_t *nions, const len_t iIon,
    const len_t Z0, const len_t rOffset
) {
    this->SetCSElements_internal(
        [&vec,&nions](const len_t i, const len_t j, const real_t v) { vec[i] += v*nions[j]; },
        iIon, Z0, rOffset
    );
}

void TritiumRateTerm::SetCSElements_internal(
    std::function<void(const len_t, const len_t, const real_t)> f,
    const len_t, const len_t, const len_t rOffset
) {
    const real_t *tritiumRate = REFluid->GetTritiumRunawayRate();

    const len_t nr = this->grid->GetNr();
    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const len_t n1 = this->grid->GetMomentumGrid(ir)->GetNp1();
        const len_t n2 = this->grid->GetMomentumGrid(ir)->GetNp2();

        const len_t xiIndex = this->GetXiIndexForEDirection(ir);
        const real_t V = this->GetVolumeScaleFactor(ir);

        // First index is into n_RE, or a distribution function...
        // Second index is into nions...
        f(offset + n1*xiIndex, rOffset+ir, scaleFactor * tritiumRate[ir] * V);

        offset += n1*n2;
    }
}

