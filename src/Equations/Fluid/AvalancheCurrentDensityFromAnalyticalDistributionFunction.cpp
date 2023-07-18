#include "DREAM/Equations/Fluid/AvalancheCurrentDensityFromAnalyticalDistributionFunction.hpp"

using namespace DREAM;

/*
 * Constructor
 */
AvalancheCurrentDensityFromAnalyticalDistributionFunction::AvalancheCurrentDensityFromAnalyticalDistributionFunction(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, RunawayFluid *rf, real_t sf
) : FVM::DiagonalComplexTerm(g,u), REFluid(rf), scaleFactor(sf) {

    id_Efield = this->unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);

    gsl_ad_w = gsl_integration_workspace_alloc(1000);
}

AvalancheCurrentDensityFromAnalyticalDistributionFunction::~AvalancheCurrentDensityFromAnalyticalDistributionFunction() {
    gsl_integration_workspace_free(gsl_ad_w);
}

void AvalancheCurrentDensityFromAnalyticalDistributionFunction::SetWeights() {
    real_t *efield = unknowns->GetUnknownData(id_Efield);
    const real_t *FSA_B = this->grid->GetRadialGrid()->GetFSA_B();
    for(len_t ir=0; ir<nr; ir++){
        real_t sgn = (efield[ir] > 0) - (efield[ir] < 0);
        real_t u_re = evaluateAvalancheRunawaysMeanVelocity(ir);
        weights[ir] = scaleFactor * sgn * Constants::ec * u_re / FSA_B[ir];
    }
}

void AvalancheCurrentDensityFromAnalyticalDistributionFunction::SetDiffWeights(len_t , len_t ) {
    // weights for Jacobian
}


/**
* Parameter struct containing integrand parameters which is passed to a GSL function.
*/
struct integrandHesslowParams {
    real_t pceff;
    real_t gammaTilde;
    real_t tauRel;
};

real_t AvalancheCurrentDensityFromAnalyticalDistributionFunction::integrandHesslow(real_t w, void *params){
    struct integrandHesslowParams *p = (struct integrandHesslowParams *)params;
    real_t mw = 1 - w;  // check if zero?
    real_t pmw = (p->pceff) + w / mw;
    real_t gt = (p->gammaTilde) * (p->tauRel);
    return gt * pmw / ( mw * mw * sqrt(1 + pmw * pmw) ) * exp( - gt * w / mw );
}

real_t AvalancheCurrentDensityFromAnalyticalDistributionFunction::evaluateAvalancheRunawaysMeanVelocity(len_t ir) {
    // some if statement...
    struct integrandHesslowParams params = {
        REFluid->GetEffectiveCriticalRunawayMomentum(ir),
        REFluid->GetAvalancheGrowthRateDividedByEMinEceff(ir),
        REFluid->GetElectronCollisionTimeRelativistic(ir)
    };

    gsl_function F;
    F.function = &(AvalancheCurrentDensityFromAnalyticalDistributionFunction::integrandHesslow);
    F.params = &params;
    real_t integral, error;

    real_t epsabs = 0, epsrel = 1e-8, lim = gsl_ad_w->limit;
    gsl_integration_qag(&F, 0, 1, epsabs, epsrel, lim, QAG_KEY, gsl_ad_w, &integral, &error);
    return integral;
}






// real_t AvalancheCurrentDensityFromAnalyticalDistributionFunction::integrandRosenbluthPutvinski(real_t s, void *params){
//     // to be implemented...
//     return 0;
// }
