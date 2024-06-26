/**
 * Implementation of the kinetic Compton scattering
 * source term, which takes the quadratic form
 *     T = S(p) * n_e(r,t)
 * where S is only a function of phase-space coordinates.
 */

#include "DREAM/Equations/Kinetic/ComptonSource.hpp"
#include "DREAM/Constants.hpp"

using namespace DREAM;

/**
 * Constructor.
 */

ComptonSource::ComptonSource(
    FVM::Grid *kineticGrid, FVM::UnknownQuantityHandler *u, FVM::Interpolator1D *comptonPhotonFlux, 
    real_t integratedComptonSpectrum, real_t C1, real_t C2, real_t C3, 
    real_t pLower, real_t scaleFactor, SourceMode sm, RunawayFluid *REFluid
) : FluidSourceTerm(kineticGrid, u), comptonPhotonFlux(comptonPhotonFlux), 
    integratedComptonSpectrum(integratedComptonSpectrum), C1(C1), C2(C2), C3(C3),
    pLower(pLower), scaleFactor(scaleFactor), sourceMode(sm), REFluid(REFluid)
{
    SetName("ComptonSource");
    this->id_ntot = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    this->id_Eterm = this->unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    this->id_ntot = this->unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    
    source = new real_t[kineticGrid->GetNCells()];
    
    this->limit = 1000;
    this->wp = gsl_integration_workspace_alloc(limit);
    this->wpOut = gsl_integration_workspace_alloc(limit);
    
    len_t offset = 0;
    this->photonFlux = this->comptonPhotonFlux->Eval(0)[0];
    
    
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++){
            for(len_t j=0; j<n2[ir]; j++){
                source[offset + j*n1[ir] + i] = EvaluateSource(ir,i,j);
            }
        }
        offset += n1[ir]*n2[ir];
    }
}

/**
 * Evaluates and sets the source for a given time (where only time-varying factor is photon flux)
 */
void ComptonSource::Rebuild(const real_t t, const real_t, FVM::UnknownQuantityHandler*) {
    this->photonFlux = this->comptonPhotonFlux->Eval(t)[0];
    this->FluidSourceTerm::Rebuild(t, 0, nullptr);
}

/**
 * Integrand of integral over photon energy for Compton Scattering, used to evaluate the number of electrons 
 * generated at each cell on the hot-tail grid. Derived from eq (23) of "Formation and termination of 
 * runaway beams in ITER disruptions" by M. Solis (2017).
 */
real_t ComptonSource::innerIntegrand(real_t Eg, void * params){
    real_t r_e = Constants::r0;
    real_t mc2 = Constants::mc2inEV;
    
    struct innerintparams* iparams = reinterpret_cast<struct innerintparams*>(params);
    real_t integratedComptonSpectrum = iparams->intConst;
    real_t C1 = iparams->c1;
    real_t C2 = iparams->c2;
    real_t C3 = iparams->c3;

    real_t z = (log(mc2 * Eg/1e6) + C1) / C2 + C3 * (mc2 * Eg/1e6) * (mc2 * Eg/1e6);
    real_t GgEg = 1 / integratedComptonSpectrum * exp(- exp(-z) - z + 1); // * photonFlux (in EvaluateSoure)
    
    real_t p = iparams->p;
    real_t g = sqrt(p*p + 1);
    real_t cost = 1 - 1 / Eg * (g - 1) / (Eg + 1 - g);
    real_t dsdO = r_e*r_e / 2 * 1 / (Eg*Eg) * (Eg / (Eg + 1 - g) + (Eg + 1 - g) / Eg - 1 + cost*cost);
    
    return GgEg * dsdO;
}

/**
 * Integrand of integral over photon energy for Compton Scattering, used to evaluate all REs generated above 
 * the hot-tail grid or above the critical momentum. 
 * In "Formation and termination of runaway beams in ITER disruptions" by M. Solis (2017), eq (23) with eq (29)
 */
real_t ComptonSource::fluidIntegrand(real_t Eg, void * params){
    real_t r_e = Constants::r0;
    real_t mc2 = Constants::mc2inEV;
    
    struct innerintparams* iparams = reinterpret_cast<struct innerintparams*>(params);
    real_t integratedComptonSpectrum = iparams->intConst;
    real_t C1 = iparams->c1;
    real_t C2 = iparams->c2;
    real_t C3 = iparams->c3;

    real_t z = (log(mc2 * Eg/1e6) + C1) / C2 + C3 * (mc2 * Eg/1e6) * (mc2 * Eg/1e6);
    real_t GgEg = 1 / integratedComptonSpectrum * exp(- exp(-z) - z + 1); // * photonFlux (in EvaluateSoure)
    
    real_t p = iparams->p;
    real_t g = sqrt(p*p + 1);
    real_t cost = 1 - 1 / Eg * (g - 1) / (Eg + 1 - g);
    real_t term1 = (Eg*Eg - 2 * Eg - 2) / (Eg*Eg*Eg) * log((1 + 2 * Eg) / (1 + Eg * (1 - cost)));
    real_t term2 = 1 / (2 * Eg) * (1 / ((1 + Eg * (1 - cost)) * (1 + Eg * (1 - cost))) - 1 / ((1 + 2 * Eg) * (1 + 2 * Eg)));
    real_t term3 = - 1 / (Eg*Eg*Eg) * (1 - Eg - (1 + 2 * Eg) / (1 + Eg * (1 - cost)) - Eg * cost);
    real_t sigma = M_PI * r_e*r_e * (term1 + term2 + term3);
    
    return GgEg * sigma;
}

/**
 * Integrand of integral to get source from Compton scattering averaged over cell on hot-tail grid
 */
real_t ComptonSource::integrand(real_t p, void * params){
    if(p == std::numeric_limits<real_t>::infinity())
        return 0.;
    real_t g = sqrt(p*p + 1);
    real_t E0 = (g-1)/2 + sqrt((g-1)*(g-1)/4+(g-1)/2);
    struct intparams* iparams = reinterpret_cast<struct intparams*>(params);
    struct innerintparams innerparams = {p, iparams->intConst, iparams->c1, iparams->c2, iparams->c3};    
    
    real_t integral;
    real_t abserr;
    gsl_function F;
    F.function = &(ComptonSource::innerIntegrand);
    F.params = &innerparams;
    gsl_integration_qagiu(&F, E0, 0., 1e-8, iparams->limit, iparams->wp, &integral, &abserr);
    
    real_t S = p / g * integral;
    return S;
}

/**
 * Evaluates the constant (only grid dependent) source-shape function S(r,p)
 */
real_t ComptonSource::EvaluateSource(len_t ir, len_t i, len_t) {
    struct intparams params = {this->limit, this->wp, this->integratedComptonSpectrum, this->C1, this->C2, this->C3};
    if(sourceMode == SOURCE_MODE_FLUID){
        return 0.;
    }
    real_t pm = grid->GetMomentumGrid(ir)->GetP1_f(i);
    real_t pp = grid->GetMomentumGrid(ir)->GetP1_f(i+1);
    if(pm == pp)
        return 0;
    real_t dp = pp-pm;
    real_t pi = (pp+pm)/2;

    if (integrand(pm, &params) < 1e-300)
        return 0.;


    real_t integral;
    real_t abserr;
    size_t neval;
    gsl_function F;
    F.function = &(ComptonSource::integrand);
    F.params = &params;
    
    gsl_integration_qng(&F, pm, pp, 0, 1e-8, &integral, &abserr, &neval);
    return scaleFactor * this->photonFlux / (2 * dp * pi*pi) * integral;
}

/**
 * Returns the source at grid point (ir,i,j).
 */
real_t ComptonSource::GetSourceFunction(len_t ir, len_t i, len_t j){
    struct intparams params = {this->limit, this->wp, this->integratedComptonSpectrum, this->C1, this->C2, this->C3};
    struct intparams paramsOut = {this->limit, this->wpOut, this->integratedComptonSpectrum, this->C1, this->C2, this->C3};
    if(sourceMode == SOURCE_MODE_FLUID){
        real_t pc = REFluid->GetEffectiveCriticalRunawayMomentum(ir);
        if (pc == std::numeric_limits<real_t>::infinity())
            return 0.;
        if (pc < pLower)
            pc = pLower;
        return scaleFactor * this->photonFlux * EvaluateTotalComptonNumber(pc, &params, &paramsOut);
    }
    
    len_t offset = 0;
    
    for(len_t iir = 0; iir < ir; iir++)
        offset += n1[iir]*n2[iir];
    len_t ind = offset + j*n1[ir] + i;
    real_t S = source[ind];
    return S;
}

/**
 * Returns the source function at (ir,i,j) differentiated with respect to the unknown x_derivId at (ir,i,j)
 */
real_t ComptonSource::GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId){
    if(derivId==id_ntot) 
        return GetSourceFunction(ir,i,j);
    else if(sourceMode == SOURCE_MODE_FLUID){
        if(derivId==id_Eterm) {
            struct intparams params = {this->limit, this->wp, this->integratedComptonSpectrum, this->C1, this->C2, this->C3};
            real_t pc = REFluid->GetEffectiveCriticalRunawayMomentum(ir);
            if (pc == std::numeric_limits<real_t>::infinity())
                return 0.;
            real_t Eeff = REFluid->GetEffectiveCriticalField(ir);
            real_t Eterm = unknowns->GetUnknownData(id_Eterm)[ir];
            real_t sgnE = (Eterm>0) - (Eterm<0);
            return -1/2 * (-integrand(pc, &params)) * pc * sgnE/( fabs(Eterm) - Eeff ) ;
        }
        else if(derivId==id_ntot) {
            struct intparams params = {this->limit, this->wp, this->integratedComptonSpectrum, this->C1, this->C2, this->C3};
            real_t pc = REFluid->GetEffectiveCriticalRunawayMomentum(ir);
            if (pc == std::numeric_limits<real_t>::infinity())
                return 0.;
            real_t ntot = unknowns->GetUnknownData(id_ntot)[ir];
            return 1/2 * (-integrand(pc, &params)) * pc / ntot;
        }
        else
            return 0;
    } 
    else
        return 0;
}

/**
 * Evaluates the number of electrons generated between pLower and pUpper, used for evaluating source outside 
 * of the hot-tail grid
 */
real_t ComptonSource::EvaluateTotalComptonNumber(real_t pLower, intparams * params, intparams * paramsOut, real_t pUpper){
    real_t integral;
    real_t abserr;
    size_t neval;
    gsl_function F;
    if(pUpper == std::numeric_limits<real_t>::infinity()) {
        struct intparams *iparamsOut = reinterpret_cast<struct intparams*>(paramsOut);
        F.function = &(fluidIntegrand);
        struct innerintparams innerparams = {pLower, iparamsOut->intConst, iparamsOut->c1, iparamsOut->c2, iparamsOut->c3};    
        F.params = &innerparams;
        real_t g = sqrt(pLower*pLower + 1);
        real_t E0 = (g-1)/2 + sqrt((g-1)*(g-1)/4+(g-1)/2);
        gsl_integration_qagiu(&F, E0, 0., 1e-8, iparamsOut->limit, iparamsOut->wp, &integral, &abserr);
        return integral;
    } else {
        F.function = &(integrand);
        F.params = params;
        gsl_integration_qng(&F, pLower, pUpper, 0, 1e-8, &integral, &abserr, &neval);
        return 2 * M_PI * integral;
    }
}
