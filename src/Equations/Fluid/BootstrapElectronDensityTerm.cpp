/**
 * Implementation of bootstrap current term driven by the electron density
 * gradient.
 */

#include "DREAM/Equations/Fluid/BootstrapElectronDensityTerm.hpp"
#include "DREAM/Constants.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
BootstrapElectronDensityTerm::BootstrapElectronDensityTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u,
    BootstrapCurrent *bs, IonHandler *ih, real_t sf
) : BootstrapEquationTerm(g, u, ih, bs, 0, sf) {

    SetUnknownID(id_ncold);

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
real_t BootstrapElectronDensityTerm::GetCoefficient(len_t ir) {

    real_t pre = bs->getConstantPrefactor(ir);
    real_t l31 = bs->getCoefficientL31(ir);

    return pre * l31 * bs->p[ir] / bs->n[ir];
}

real_t BootstrapElectronDensityTerm::GetPartialCoefficient(len_t ir, len_t derivId, len_t jzs, len_t /* jZ */) {

    real_t pre = bs->getConstantPrefactor(ir);
    real_t l31 = bs->getCoefficientL31(ir);
    real_t dl31 = bs->evaluatePartialCoefficientL31(ir, derivId, jzs);

    real_t dCoefficient = dl31;
    if (!bs->includeIonTemperatures) { 
        dCoefficient *= bs->Tcold[ir];
        if (derivId == id_Tcold)
            dCoefficient += l31;
    } else {
        dCoefficient *= bs->p[ir] / bs->n[ir];
        if (derivId == id_ncold)
            dCoefficient += l31 * (bs->Tcold[ir] - bs->p[ir] / bs->n[ir]) / bs->n[ir];
        else if (derivId == id_Tcold)
            dCoefficient += l31 * bs->ncold[ir] / bs->n[ir];
        else if (derivId == id_Ni)
            dCoefficient -= l31 * bs->p[ir] / (bs->n[ir] * bs->n[ir]);
        else if (derivId == id_Wi)
                dCoefficient += l31 / (1.5 * bs->n[ir] * Constants::ec);
    }
    return pre * dCoefficient;
}
