#include "DREAM/Equations/Fluid/RunawayCurrentDensityFromAnalyticalDistributionFunction.hpp"
#include <iostream>
using namespace DREAM;

/**
 * Constructor.
 */
RunawayCurrentDensityFromAnalyticalDistributionFunction::RunawayCurrentDensityFromAnalyticalDistributionFunction(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, RunawayFluid *rf, real_t sf
) : FVM::DiagonalComplexTerm(g,u), REFluid(rf), scaleFactor(sf) {

    id_Efield = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    id_ntot   = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);

    gsl_w = gsl_integration_workspace_alloc(GSL_WORKSPACE_SIZE);

    AddUnknownForJacobian(unknowns, id_Efield);
    AddUnknownForJacobian(unknowns, id_ntot);
}

/**
 * Destructor.
 */
RunawayCurrentDensityFromAnalyticalDistributionFunction::~RunawayCurrentDensityFromAnalyticalDistributionFunction() {
    gsl_integration_workspace_free(gsl_w);
    if (dPCrit != nullptr)
        delete [] dPCrit;
}

/**
 * Set the weights of this term.
 */
void RunawayCurrentDensityFromAnalyticalDistributionFunction::SetWeights() {
    Efield = unknowns->GetUnknownData(id_Efield);
    for(len_t ir=0; ir < nr; ir++){
        real_t sgn = (Efield[ir] > 0) - (Efield[ir] < 0);
        real_t u_re = evaluateMeanSpeed(ir);
        const real_t FSA_B = this->grid->GetRadialGrid()->GetFSA_B(ir);
        weights[ir] = scaleFactor * sgn * constPreFactor * u_re / FSA_B;
    }
}

/**
 * Set the weights for the weight Jacobian matrix of this term.
 *
 * Assumes that the only dependence of u_re to any unknown quantities is in the lower
 * integral boundary p = pCrit is captured via the approximation pCrit ~ 1/sqrt{E - Eceff},
 * with Eceff ~ ntot.
 */
void RunawayCurrentDensityFromAnalyticalDistributionFunction::SetDiffWeights(len_t derivId, len_t /*nMultiples*/ ) {
    if ( !((derivId == id_Efield) || (derivId == id_ntot)) ) {
        for (len_t ir = 0; ir < nr; ir++)
            diffWeights[ir] = 0;
    } else {
        if (dPCrit == nullptr)
            dPCrit = new real_t[nr];

        REFluid->evaluatePartialContributionCriticalREMomentum(dPCrit, derivId);
        struct integrandParams params;
        for (len_t ir = 0; ir < nr; ir++) {
            if (dPCrit[ir] == 0)
                diffWeights[ir] = 0;
            else {
                params = {ir, Efield[ir], REFluid};
                const real_t FSA_B = grid->GetRadialGrid()->GetFSA_B(ir);
                diffWeights[ir] = - scaleFactor * constPreFactor * integrand(0, &params) * dPCrit[ir] / FSA_B;  // Leibniz integral rule
            }
        }
    }
}

/**
 * Returns the integrand appearing in the evaluation of the mean RE speed.
 */
real_t RunawayCurrentDensityFromAnalyticalDistributionFunction::integrand(real_t w, void *params) {
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
real_t RunawayCurrentDensityFromAnalyticalDistributionFunction::evaluateMeanSpeed(len_t ir) {

    if (std::isinf(REFluid->GetEffectiveCriticalRunawayMomentum(ir)))
        return 1;

    gsl_function gsl_func;
    gsl_func.function = &(RunawayCurrentDensityFromAnalyticalDistributionFunction::integrand);

    struct integrandParams params = {ir, Efield[ir], REFluid};
    gsl_func.params = &params;

    real_t integral, error;
    real_t epsabs = 0, epsrel = 1e-8, lim = gsl_w->limit;
    gsl_integration_qag(&gsl_func, 0, 1, epsabs, epsrel, lim, QAG_KEY, gsl_w, &integral, &error);

    return integral;
}
