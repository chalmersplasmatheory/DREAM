#include "DREAM/Equations/Fluid/AvalancheCurrentDensityFromAnalyticalDistributionFunction.hpp"
#include <iostream>
using namespace DREAM;

/**
 * Constructor.
 */
AvalancheCurrentDensityFromAnalyticalDistributionFunction::AvalancheCurrentDensityFromAnalyticalDistributionFunction(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, RunawayFluid *rf, real_t sf
) : FVM::DiagonalComplexTerm(g,u), REFluid(rf), scaleFactor(sf) {
    id_Efield = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    gsl_w = gsl_integration_workspace_alloc(GSL_WORKSPACE_SIZE);
}

/**
 * Destructor.
 */
AvalancheCurrentDensityFromAnalyticalDistributionFunction::~AvalancheCurrentDensityFromAnalyticalDistributionFunction() {
    gsl_integration_workspace_free(gsl_w);
}

/**
 * Set the weights of this term.
 */
void AvalancheCurrentDensityFromAnalyticalDistributionFunction::SetWeights() {
    Efield = unknowns->GetUnknownData(id_Efield);
    const real_t *FSA_B = this->grid->GetRadialGrid()->GetFSA_B();
    for(len_t ir=0; ir<nr; ir++){
        real_t sgn = (Efield[ir] > 0) - (Efield[ir] < 0);
        real_t u_re = evaluateMeanSpeed(ir);
        weights[ir] = scaleFactor * sgn * Constants::c * Constants::ec * u_re / FSA_B[ir];
    }
}

/**
 * Set the weights for the Jacobian matrix of this term.
 */
void AvalancheCurrentDensityFromAnalyticalDistributionFunction::SetDiffWeights(len_t , len_t ) {
    // weights for Jacobian
}


/**
* Parameter struct used for the evaluation of mean RE speed integral.
*/
struct integrandParams {
    len_t ir;
    real_t Efield;
    RunawayFluid *REFluid;
};

/**
 * Returns the integrand appearing in the evaluation of the mean RE speed.
 */
real_t AvalancheCurrentDensityFromAnalyticalDistributionFunction::integrand(real_t w, void *params) {
    struct integrandParams *p = (struct integrandParams *)params;

    RunawayFluid *rf = p->REFluid;
    real_t EMinusEceff = (fabs(p->Efield) - rf->GetEffectiveCriticalField(p->ir)) * Constants::ec / (Constants::me * Constants::c);
    real_t beta = rf->GetAvalancheGrowthRate(p->ir) / EMinusEceff;

    real_t mw = 1 - w;  // check if zero?
    real_t pmw = rf->GetEffectiveCriticalRunawayMomentum(p->ir) + w / mw;
    return beta * pmw / ( mw * mw * sqrt(1 + pmw * pmw) ) * exp( - beta * w / mw );
}

/**
 * Calculates the mean RE speed assuming an analytical RE distribution function, based on Eq. (4.2) in Svensson et al. (JPP 2020).
 */
real_t AvalancheCurrentDensityFromAnalyticalDistributionFunction::evaluateMeanSpeed(len_t ir) {

    if (std::isinf(REFluid->GetEffectiveCriticalRunawayMomentum(ir)))
        return 1;

    gsl_function gsl_func;
    gsl_func.function = &(AvalancheCurrentDensityFromAnalyticalDistributionFunction::integrand);

    struct integrandParams params = {ir, Efield[ir], REFluid};
    gsl_func.params = &params;

    real_t integral, error;
    real_t epsabs = 0, epsrel = 1e-8, lim = gsl_w->limit;
    gsl_integration_qag(&gsl_func, 0, 1, epsabs, epsrel, lim, QAG_KEY, gsl_w, &integral, &error);

    return integral;
}
