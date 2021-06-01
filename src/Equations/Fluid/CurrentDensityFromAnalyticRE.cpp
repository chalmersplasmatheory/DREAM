/**
 * Implementation of equation term representing the
 * current density carried by the analytic RE distribution
 */

#include "DREAM/Equations/Fluid/CurrentDensityFromAnalyticRE.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
CurrentDensityFromAnalyticRE::CurrentDensityFromAnalyticRE(           
    FVM::Grid *grid, FVM::UnknownQuantityHandler *u,
    AnalyticDistributionRE *dist, real_t sf
) : FVM::EquationTerm(grid), distRE(dist), unknowns(u),
    id_Efield(u->GetUnknownID(OptionConstants::UQTY_E_FIELD)),
    id_ncold(u->GetUnknownID(OptionConstants::UQTY_N_COLD)),
    scaleFactor(sf)
{
    AddUnknownForJacobian(unknowns,id_Efield);
    //AddUnknownForJacobian(unknowns,id_ncold);

    real_t preFactor = scaleFactor*Constants::ec*Constants::c;
    FVM::RadialGrid *rGrid = grid->GetRadialGrid();
    AveragedXiTerm = new REPitchDistributionAveragedBACoeff(
        rGrid, distRE, &(FVM::RadialGrid::BA_FUNC_XI_OVER_B), nullptr, 
        FVM::RadialGrid::BA_PARAM_XI_OVER_B, [](real_t xi0){return xi0;},
        [this,preFactor](len_t /*ir*/, real_t p){
            return preFactor * p / sqrt(1+p*p);
    });

    gsl_ad = gsl_integration_workspace_alloc(1000);
    gsl_func.function = &(currentDensityIntegrand);
    gsl_func.params = &gsl_params;

    gsl_func_partial.function = &(partialCurrentDensityIntegrand);
    gsl_func_partial.params   = &gsl_params;


    gsl_params.distRE = distRE;
    gsl_params.AveragedXiTerm = AveragedXiTerm;

    GridRebuilt();
}

/**
 * Destructor
 */
CurrentDensityFromAnalyticRE::~CurrentDensityFromAnalyticRE(){
    Deallocate();
    gsl_integration_workspace_free(gsl_ad);
    delete [] AveragedXiTerm;
}           

void CurrentDensityFromAnalyticRE::Deallocate(){
    if(currentDensity != nullptr){
        delete [] currentDensity;
        delete [] partialCurrentDensity;
    }   
}

/**
 * Called after the grid is rebuilt; (re)allocates memory
 * for all quantities
 */
bool CurrentDensityFromAnalyticRE::GridRebuilt(){
    Deallocate();
    currentDensity = new real_t[nr];
    partialCurrentDensity = new real_t[nr];
    AveragedXiTerm->GridRebuilt();
    return true;
}


real_t CurrentDensityFromAnalyticRE::partialCurrentDensityIntegrand(real_t p, void *par){
    integrandParameters *params = (integrandParameters*) par;
    len_t ir = params->ir;
    len_t derivId = params->derivId;
    len_t nMultiple = params->nMultiple;
    real_t pitchDistAverage;
    real_t partialPitchDistAverage = params->AveragedXiTerm->EvaluatePartialREPitchDistAverage(ir, p, derivId, nMultiple, &pitchDistAverage);
    real_t energyDist;
    real_t partialEnergyDist = params->distRE->evaluatePartialEnergyDistribution(ir, p, derivId, nMultiple, &energyDist);
    return p*p*(partialPitchDistAverage * energyDist + pitchDistAverage * partialEnergyDist);
}

real_t CurrentDensityFromAnalyticRE::currentDensityIntegrand(real_t p, void *par){
    integrandParameters *params = (integrandParameters*) par;
    len_t ir = params->ir;
    return p*p * params->AveragedXiTerm->EvaluateREPitchDistAverage(ir, p) * params->distRE->evaluateEnergyDistribution(ir, p);
}

/**
 * Rebuilds quantities used by this equation term.
 */
void CurrentDensityFromAnalyticRE::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) {
    real_t 
        epsabs = 0,
        epsrel = 5e-4,
        limit  = gsl_ad->limit,
        error  = 0;
    for(len_t ir=0; ir<nr; ir++){
        gsl_params.ir = ir;
        // adaptive integration over momentum from p=0 to infinity 
        gsl_integration_qagiu(&gsl_func, 0.0, epsabs, epsrel, limit, gsl_ad, &currentDensity[ir], &error);
    }
}

void CurrentDensityFromAnalyticRE::setPartialCurrentDensity(len_t derivId, len_t nMultiple){
    real_t 
        epsabs = 0,
        epsrel = 5e-4,
        limit  = gsl_ad->limit,
        error  = 0;
    gsl_params.derivId = derivId;
    gsl_params.nMultiple = nMultiple;
    for(len_t ir=0; ir<nr; ir++){
        gsl_params.ir = ir;
        // adaptive integration over momentum from p=0 to infinity 
        gsl_integration_qagiu(&gsl_func_partial, 0.0, epsabs, epsrel, limit, gsl_ad, &partialCurrentDensity[ir], &error);
    }

}

void CurrentDensityFromAnalyticRE::SetVectorElements(real_t *vec, const real_t *n_re){
    for(len_t ir=0; ir<nr; ir++)
        vec[ir] += n_re[ir]*currentDensity[ir];
}



bool CurrentDensityFromAnalyticRE::SetJacobianBlock(
    len_t uqtyId, len_t derivId, FVM::Matrix *jac, const real_t *n_re
) {
    bool contrib = false;
    if(uqtyId==derivId){
        contrib = true;
        SetMatrixElements(jac,nullptr);
    }

    len_t nMultiples;
    if(!HasJacobianContribution(derivId, &nMultiples))
        return contrib;
    
    for(len_t n=0; n<nMultiples; n++){
        setPartialCurrentDensity(derivId,n);
        for(len_t ir=0; ir<nr; ir++)
            jac->SetElement(ir,nr*n + ir, n_re[ir]*partialCurrentDensity[ir]);
    }

    return true;
}

void CurrentDensityFromAnalyticRE::SetMatrixElements(FVM::Matrix *mat, real_t*) {
    for(len_t ir=0; ir<nr; ir++){
        mat->SetElement(ir,ir,currentDensity[ir]);
    }
}

