/**
 * Implementation of bootstrap current term driven by the electron density
 * gradient.
 */
#include "DREAM/Equations/Fluid/BootstrapElectronDensityTerm.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
BootstrapElectronDensityTerm::BootstrapElectronDensityTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, BootstrapCurrent *bs, IonHandler *ih, real_t sf
) : BootstrapEquationTerm(g, u, ih, sf), bs(bs) {

    SetUnknown(id_ncold);

    SetName("BootstrapElectronDensityTerm");
    AddUnknownForJacobian(u, id_ncold);
    AddUnknownForJacobian(u, id_Tcold);
    AddUnknownForJacobian(u, id_Ni);
    AddUnknownForJacobian(u, id_ions);
    if (bs->includeIonTemperatures)
        AddUnknownForJacobian(u, id_Wi);
}

/**
 * Coefficient = A * L31 * p / n
 */
real_t BootstrapElectronDensityTerm::GetCoefficient(len_t ir, len_t /* iZ */) {
    real_t coefficient;
    if (bs->includeIonTemperatures)
        coefficient = bs->p[ir] / bs->n[ir];
    else
        coefficient = bs->Tcold[ir];
    return bs->constantPrefactor[ir] * bs->coefficientL31[ir] * coefficient;
}

real_t BootstrapElectronDensityTerm::GetPartialCoefficient(
    len_t ir, len_t derivId, len_t index, len_t /* iZ */
) {
    real_t dCoefficient = evaluatePartialCoefficientL31(ir, derivId, index)
    if (!bs->includeIonTemperatures) {
        dCoefficient *= bs->Tcold[ir];
        dCoefficient += bs->coefficientL31[ir] * (real_t)(derivId == id_Tcold);
    } else {
        dCoefficient *= bs->p[ir] / bs->n[ir];
        switch (derivId) {
            case id_ncold:
                dCoefficient += bs->coefficientL31[ir] * (bs->Tcold[ir] - bs->p[ir] / bs->n[ir]) / bs->n[ir];
                break;
            case id_Tcold:
                dCoefficient += bs->coefficientL31[ir] * bs->ncold[ir] / bs->n[ir];
                break;
            case id_Ni:
                dCoefficient -= bs->coefficientL31[ir] / (bs->n[ir] * bs->n[ir]);
                break;
            case id_Wi:
                dCoefficient += bs->coefficientL31[ir] * 2. / (3. * bs->n[ir]);
        }
    }
    return bs->constantPrefactor[ir] * dCoefficient;
}
