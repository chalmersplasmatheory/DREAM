/**
 * Implementation of bootstrap current term driven by the ion thermal energy
 * gradients.
 */
#include "DREAM/Equations/Fluid/BootstrapIonThermalEnergyTerm.hpp"
#include "DREAM/Constants.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
BootstrapIonThermalEnergyTerm::BootstrapIonThermalEnergyTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u,
    BootstrapCurrent *bs, IonHandler *ih, real_t sf
) : BootstrapEquationTerm(g, u, ih, bs, sf) {

    SetUnknownID(id_Wi);

    SetName("BootstrapIonThermalEnergyTerm");
    AddUnknownForJacobian(u, id_ncold);
    AddUnknownForJacobian(u, id_Tcold);
    AddUnknownForJacobian(u, id_Ni);
    AddUnknownForJacobian(u, id_ions);
}

/**
 * If ion temperature evolution is included, then for each ion:
 *      Coefficient = A * L31 * (1 + alpha) * 2/3
 *
 * If not included, ie. Ti = Tcold, then this term is not used!
 *
 * Note, we divide with the electron charge because W_i is given in SI units!
 */
real_t BootstrapIonThermalEnergyTerm::GetCoefficient(len_t ir, len_t /* iZ */) {

    real_t pre = bs->getConstantPrefactor(ir);
    real_t l31 = bs->getCoefficientL31(ir);
    real_t alpha = bs->getCoefficientAlpha(ir);

    return pre * l31 * 2./3. * ( 1. + alpha ) / Constants::ec;
}

/**
 * Partial derivative of Coefficient.
 */
real_t BootstrapIonThermalEnergyTerm::GetPartialCoefficient(len_t ir, len_t derivId, len_t index, len_t iZ) {

    real_t pre = bs->getConstantPrefactor(ir);
    real_t l31 = bs->getCoefficientL31(ir);
    real_t alpha = bs->getCoefficientAlpha(ir);
    real_t dl31 = bs->evaluatePartialCoefficientL31(ir, derivId, index);
    real_t dalpha = bs->evaluatePartialCoefficientAlpha(ir, derivId, index, iZ);

    return pre * 2./3. *  ( dl31 * ( 1. + alpha ) + l31 * dalpha );
}
