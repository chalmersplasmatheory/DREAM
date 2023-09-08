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
    id_ncold  = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);

    FSA_B = this->grid->GetRadialGrid()->GetFSA_B();

    u_re = new real_t[nr];
    beta = new real_t[nr];
    gsl_w = gsl_integration_workspace_alloc(GSL_WORKSPACE_SIZE);

    AddUnknownForJacobian(unknowns, id_Efield);
    AddUnknownForJacobian(unknowns, id_ntot);
    AddUnknownForJacobian(unknowns, id_ncold);
}

/**
 * Destructor.
 */
RunawayCurrentDensityFromAnalyticalDistributionFunction::~RunawayCurrentDensityFromAnalyticalDistributionFunction() {
    gsl_integration_workspace_free(gsl_w);
    delete [] u_re;
    delete [] beta;
    if (dPCrit != nullptr)
        delete [] dPCrit;
}

/**
 * Set the weights of this term.
 */
void RunawayCurrentDensityFromAnalyticalDistributionFunction::SetWeights() {
    evaluateMeanSpeed();
    for(len_t ir = 0; ir < nr; ir++){
        real_t sgn = (Efield[ir] > 0) - (Efield[ir] < 0);
        weights[ir] = scaleFactor * sgn * constPreFactor * u_re[ir] / FSA_B[ir];
    }
}

/**
 * Set the weights for the weight Jacobian matrix of this term. d(j_re) / d(uqty)
 *
 * Assumes that the only dependence of u_re to any unknown quantities is in the lower
 * integral boundary p = pCrit is captured via the approximation pCrit ~ sqrt{ntot/(E-Eceff)},
 * with Eceff ~ ntot.
 */
void RunawayCurrentDensityFromAnalyticalDistributionFunction::SetDiffWeights(len_t derivId, len_t /*nMultiples*/ ) {

    if ( !((derivId == id_Efield) || (derivId == id_ntot) || (derivId == id_ncold)) ) {
        for (len_t ir = 0; ir < nr; ir++)
            diffWeights[ir] = 0;
    }
    if (derivId == id_ncold) {
        real_t *n_cold = unknowns->GetUnknownData(id_ncold);
        for (len_t ir = 0; ir < nr; ir++) {
            real_t pCrit = REFluid->GetEffectiveCriticalRunawayMomentum(ir);
            if (isinf(pCrit))
                diffWeights[ir] = 0;
            else {
                real_t du_re = - u_re[ir] / n_cold[ir] * ( 1 + beta[ir] * pCrit );
                diffWeights[ir] = scaleFactor * constPreFactor * du_re / FSA_B[ir];
            }
        }
    } else {
        if (dPCrit == nullptr)
            dPCrit = new real_t[nr];

        REFluid->evaluatePartialContributionCriticalREMomentum(dPCrit, derivId);
        const real_t *pCrit = REFluid->GetEffectiveCriticalRunawayMomentum();
        struct integrandParams params;

        if (derivId == id_Efield) {
            for (len_t ir = 0; ir < nr; ir++) {
                if (dPCrit[ir] == 0)
                    diffWeights[ir] = 0;
                else {
                    params = {beta[ir], pCrit[ir]};
                    real_t du_re = - dPCrit[ir] * ( integrand(0, &params) - beta[ir] * u_re[ir] );
                    diffWeights[ir] = scaleFactor * constPreFactor * du_re / FSA_B[ir];
                }
            }
        } else { /* id_ntot */
            real_t *n_tot = unknowns->GetUnknownData(id_ntot);
            for (len_t ir = 0; ir < nr; ir++) {
                if (dPCrit[ir] == 0)
                    diffWeights[ir] = 0;
                else {
                    params = {beta[ir], pCrit[ir]};
                    real_t du_re = - dPCrit[ir] * ( integrand(0, &params) - 3 * u_re[ir] * pCrit[ir] * beta[ir] ) + u_re[ir] / n_tot[ir];
                    diffWeights[ir] = scaleFactor * constPreFactor * du_re / FSA_B[ir];
                }
            }
        }
    }
}

/**
 * Returns the integrand appearing in the evaluation of the mean RE speed.
 */
real_t RunawayCurrentDensityFromAnalyticalDistributionFunction::integrand(real_t w, void *params) {
    struct integrandParams *p = (struct integrandParams *)params;
    real_t mw = 1 - w;  // check if zero?
    real_t pmw = p->pCrit + w / mw;
    return p->beta * pmw / ( mw * mw * sqrt(1 + pmw * pmw) ) * exp( - p->beta * w / mw );
}

/**
 * Calculates the mean RE speed assuming an analytical RE distribution function, based on Eq. (4.2) in Svensson et al. (JPP 2020).
 */
void RunawayCurrentDensityFromAnalyticalDistributionFunction::evaluateMeanSpeed() {

    Efield = unknowns->GetUnknownData(id_Efield);

    const real_t *pCrit = REFluid->GetEffectiveCriticalRunawayMomentum();
    const real_t *Eceff = REFluid->GetEffectiveCriticalField();
    const real_t *gamma = REFluid->GetAvalancheGrowthRate();

    gsl_function gsl_func;

    for (len_t ir = 0; ir < nr; ir++) {

        if (std::isinf(pCrit[ir]))
            u_re[ir] = 1;   // is this the correct option?
        else {

            // beta is used both here (call from SetWeights) and in SetDiffWeights.
            real_t EMinusEceff = ( fabs(Efield[ir]) - Eceff[ir] ) * Constants::ec / ( Constants::me * Constants::c );
            beta[ir] = gamma[ir] / EMinusEceff;

            struct integrandParams params = {beta[ir], pCrit[ir]};
            gsl_func.params = &params;
            gsl_func.function = &(RunawayCurrentDensityFromAnalyticalDistributionFunction::integrand);

            real_t error;
            real_t epsabs = 0, epsrel = 1e-8, lim = gsl_w->limit;
            gsl_integration_qag(&gsl_func, 0, 1, epsabs, epsrel, lim, QAG_KEY, gsl_w, &u_re[ir], &error);
        }
    }
}
