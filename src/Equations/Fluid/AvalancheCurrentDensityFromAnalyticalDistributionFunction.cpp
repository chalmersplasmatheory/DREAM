#include "DREAM/Equations/AvalancheCurrentDensityFromAnalyticalDistributionFunction.hpp"

using namespace DREAM;

/*
 * Constructor
 */
AvalancheCurrentDensityFromAnalyticalDistributionFunction::AvalancheCurrentDensityFromAnalyticalDistributionFunction(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, RunawayFluid *rf, real_t sf
) : FVM::DiagonalComplexTerm(g,u), REFluid(rf), scaleFactor(sf) {


    id_ncold  = this->unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_ntot   = this->unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    id_Efield = this->unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);


    ncold  = new real_t[nr];
    ntot   = new real_t[nr];
    efield = new real_t[nr];
    pceff  = new real_t[nr];
    nuDnuS = new real_t[nr];



    gsl_ad_w = gsl_integration_workspace_alloc(1000);
}

AvalancheCurrentDensityFromAnalyticalDistributionFunction::~AvalancheCurrentDensityFromAnalyticalDistributionFunction() {

    delete [] ncold;
    delete [] ntot;
    delete [] efield;
    delete [] pceff;
    delete [] nuDnuS;

    gsl_integration_workspace_free(gsl_ad_w);
}

void AvalancheCurrentDensityFromAnalyticalDistributionFunction::SetWeights() {
    // weights for matrix

    ncold  = unknowns->GetUnknownData(id_ncold);
    ntot   = unknowns->GetUnknownData(id_ntot);
    efield = unknowns->GetUnknownData(id_Efield);

    const real_t *FSA_B = this->grid->GetRadialGrid()->GetFSA_B();
    for(len_t ir=0; ir<nr; ir++){
        real_t sgn = (Efield[ir] > 0) - (Efield[ir] < 0);
        weights[ir] = sgn * Constants::ec * Constants::c / FSA_B[ir];
    }
}

void AvalancheCurrentDensityFromAnalyticalDistributionFunction::SetDiffWeights(len_t derivId, len_t nMultiples) {
    // weights for Jacobian
}



real_t AvalancheCurrentDensityFromAnalyticalDistributionFunction::integrandHesslow(real_t w, void *params){
    struct integrandHesslowParams *p = (struct integrandHesslowParams *)params;
    real_t factor = (p->ncold) / ((p->ntot) * (p->lnLambda));
    real_t pceff = (p->pceff);
    real_t collfreqs = 1 / sqrt( 4 + (p->nuDnuS) );
    real_t mw = 1 - w;  // check if zero?
    real_t ex = exp( -factor * collfreqs * w / mw );
    real_t pmw = pceff + w / mw;
    return factor * collfreqs * pmw * ex / ( mw * mw * sqrt( 1 + pmw * pmw ));
}

real_t AvalancheCurrentDensityFromAnalyticalDistributionFunction::evaluateAvalancheRunawaysMeanVelocity(len_t ir) {
    // some if statement...
    struct integrandHesslowParams params = {ncold[ir], ntot[ir], lnLambda[ir], pceff[ir], nuDnuS[ir]};

    gsl_function F;
    F.function = &(AvalancheCurrentDensityFromAnalyticalDistributionFunction::integrandHesslow);
    F.params = &params;
    real_t integral, error;

    real_t epsabs = 0, epsrel = 1e-8, lim = gsl_ad_w->limit;
    gsl_integration_qag(&F, 0, 1, epsabs, epsrel, lim, QAG_KEY, gsl_ad_w, &integral, &error);
    return integral;
}






real_t AvalancheCurrentDensityFromAnalyticalDistributionFunction::integrandRosenbluthPutvinski(real_t s, void *params){
    // to be implemented...
    return 0;
}
