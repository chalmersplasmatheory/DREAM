/**
 * Implementation of bootstrap current term driven by the electron temperature
 * gradient.
 */
#include "DREAM/Equations/Fluid/BootstrapElectronTemperatureTerm.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
BootstrapElectronTemperatureTerm::BootstrapElectronTemperatureTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, BootstrapCurrent *bs,
    IonHandler *ih, OptionConstants::eqterm_bootstrap_bc bc, real_t sf
) : BootstrapEquationTerm(g, u, ih, bs, bc, sf) {

    SetUnknownID(id_Tcold);

    SetName("BootstrapElectronTemperatureTerm");
    AddUnknownForJacobian(u, id_ncold);
    AddUnknownForJacobian(u, id_Tcold);
    AddUnknownForJacobian(u, id_Ni);
    AddUnknownForJacobian(u, id_ions);
}

/**
 * If ion temperature evolution is included:
 *      Coefficient = A * (L31 + L32) * ncold
 *
 * If not included, ie. Ti = Tcold, then:
 *      Coefficient = A * [ (L31 + L32) * ncold + L31 * (1 + alpha) * ionsum(ni) ]
 */
real_t BootstrapElectronTemperatureTerm::GetCoefficient(len_t ir, len_t /* iZ */) {
    real_t coefficient = ( bs->coefficientL31[ir] + bs->coefficientL32[ir] ) * bs->ncold[ir];
    if (!bs->includeIonTemperatures) {
        real_t nitot = 0;
        for (len_t i = ir; i < nr * nZ; i += nr)
            nitot += bs->Ni[i];
        coefficient += bs->coefficientL31[ir] * ( 1. + bs->coefficientAlpha[ir] ) * nitot;
    }
    return bs->constantPrefactor[ir] * coefficient;
}

/**
 * Partial derivative of Coefficient.
 *
 * ir:      radial cell grid point.
 * derivId: unknown quantity ID to differentiate with respect to.
 * iz:      ion charge state index
 * iZ:      ion species index
 */
real_t BootstrapElectronTemperatureTerm::GetPartialCoefficient(
    len_t ir, len_t derivId, len_t iz, len_t iZ
) {
    real_t dl31 = bs->evaluatePartialCoefficientL31(ir, derivId, iz);
    real_t dl32 = bs->evaluatePartialCoefficientL32(ir, derivId, iz);
    real_t dCoefficient = (dl31 + dl32) * bs->ncold[ir];
    if (derivId == id_ncold)
        dCoefficient += bs->coefficientL31[ir] + bs->coefficientL32[ir];
    if (!bs->includeIonTemperatures) {
        real_t nitot = 0;
        for (len_t i = ir; i < nr * nZ; i += nr)
            nitot += bs->Ni[i];
        dCoefficient += dl31 * ( 1. + bs->coefficientAlpha[ir] ) * nitot;

        real_t dAlpha = bs->evaluatePartialCoefficientAlpha(ir, derivId, iz, iZ);
        dCoefficient += bs->coefficientL31[ir] * dAlpha * nitot;

        if (derivId == id_Ni)
            dCoefficient += bs->coefficientL31[ir] * ( 1. + bs->coefficientAlpha[ir] );
    }
    return bs->constantPrefactor[ir] * dCoefficient;
}
