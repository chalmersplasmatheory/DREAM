#include "DREAM/Equations/Fluid/AvalancheCurrentDensityFromAnalyticalDistributionFunction.hpp"
#include <iostream>
using namespace DREAM;

/*
 * Constructor
 */
AvalancheCurrentDensityFromAnalyticalDistributionFunction::AvalancheCurrentDensityFromAnalyticalDistributionFunction(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, RunawayFluid *rf,/* OptionConstants::eqterm_fluid_runaway_current_mode *fm,*/ real_t sf
) : FVM::DiagonalComplexTerm(g,u), REFluid(rf), /*fm(fm),*/ scaleFactor(sf) {

    id_Efield = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);

    // old_u_re = new real_t[nr];

    gsl_w0 = gsl_integration_workspace_alloc(GSL_WORKSPACE_SIZE);

    // for 2d integrals, use double integration
    // if (fm == OptionConstants::EQTERM_FLUID_RUNAWAY_CURRENT_MODE_ROSENBLUTH_PUTVINSKI_MOMENT)
        // gsl_w1 = gsl_integration_workspace_alloc(GSL_WORKSPACE_SIZE);
}

AvalancheCurrentDensityFromAnalyticalDistributionFunction::~AvalancheCurrentDensityFromAnalyticalDistributionFunction() {

    // delete [] old_u_re;

    gsl_integration_workspace_free(gsl_w0);
    if (gsl_w1 != nullptr)
        gsl_integration_workspace_free(gsl_w1);
}

void AvalancheCurrentDensityFromAnalyticalDistributionFunction::SetWeights() {
    real_t *efield = unknowns->GetUnknownData(id_Efield);
    // real_t *eceff  = REFluid->GetEffectiveCriticalField();
    const real_t *FSA_B = this->grid->GetRadialGrid()->GetFSA_B();

    for(len_t ir=0; ir<nr; ir++){
        real_t sgn = (efield[ir] > 0) - (efield[ir] < 0);

        real_t u_re = evaluateAvalancheRunawaysMeanVelocity(ir);
        // if (fabs(efield[ir]) > eceff[ir])
        //     u_re *= evaluateAvalancheRunawaysMeanVelocity(ir);
        std::cout << u_re << std::endl;
        weights[ir] = scaleFactor * sgn * Constants::c * Constants::ec * u_re / FSA_B[ir];
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

real_t AvalancheCurrentDensityFromAnalyticalDistributionFunction::integrandHesslow(real_t w, void *params) {
    struct integrandHesslowParams *p = (struct integrandHesslowParams *)params;
    real_t mw = 1 - w;  // check if zero?
    real_t pmw = (p->pceff) + w / mw;
    real_t gt = (p->gammaTilde) * (p->tauRel);
    return gt * pmw / ( mw * mw * sqrt(1 + pmw * pmw) ) * exp( - gt * w / mw );
    // return gt / ( mw * mw ) * exp( - gt * w / mw );
}


// struct integrandRosenbluthPutvinskiParams {
//     real_t efield;
//     real_t ecrit;
//     real_t zeff;
// };

// real_t AvalancheCurrentDensityFromAnalyticalDistributionFunction::integrandRosenbluthPutvinski(real_t w, void *params) {
//     struct integrandRosenbluthPutvinskiParams *p = (struct integrandRosenbluthPutvinskiParams *)params;
//
//     return 0;
// }
//
// struct integrandRosenbluthPutvinski_innerParams {
//     real_t w;
//     real_t efield;
//     real_t ecrit;
//     real_t zeff;
// };


// real_t AvalancheCurrentDensityFromAnalyticalDistributionFunction::integrandRosenbluthPutvinski_inner(real_t z, void *params) {
//     return 0;
// }


real_t AvalancheCurrentDensityFromAnalyticalDistributionFunction::evaluateAvalancheRunawaysMeanVelocity(len_t ir) {
    // some if statement...

    // struct integrandRosenbluthPutvinskiParams params = {
    //     nullptr,
    //     // unknowns->GetUnknownData(id_Efield);
    //     REFluid->GetConnorHastieField_COMPLETESCREENING(ir),
    // };

    // hesslow

    real_t pceff = REFluid->GetEffectiveCriticalRunawayMomentum(ir);
    if (std::isinf(pceff)) {
        return 1;
    }

    struct integrandHesslowParams params = {
        pceff,
        REFluid->GetAvalancheGrowthRateDividedByEMinEceff(ir),
        REFluid->GetElectronCollisionTimeRelativistic(ir)
    };

    gsl_function F;
    F.function = &(AvalancheCurrentDensityFromAnalyticalDistributionFunction::integrandHesslow);
    F.params = &params;
    real_t integral, error;

    real_t epsabs = 0, epsrel = 1e-8, lim = gsl_w0->limit;
    gsl_integration_qag(&F, 0, 1, epsabs, epsrel, lim, QAG_KEY, gsl_w0, &integral, &error);
    return integral;
}
