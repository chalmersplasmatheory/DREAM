/**
 * Implementation of the kinetic Compton scattering
 * source term, which takes the quadratic form
 *     T = S(p) * n_e(r,t)
 * where S is only a function of phase-space coordinates.
 */

#include "DREAM/Equations/Kinetic/ComptonSource.hpp"
#include "DREAM/Constants.hpp"
//#include <gsl/gsl_integration.h> // ??

using namespace DREAM;

/**
 * Constructor.
 */
ComptonSource::ComptonSource(
    FVM::Grid *kineticGrid, FVM::UnknownQuantityHandler *u, real_t photonFlux, real_t pLower, real_t scaleFactor, SourceMode sm, RunawayFluid *REFluid
) : FluidSourceTerm(kineticGrid, u), photonFlux(photonFlux), pLower(pLower), scaleFactor(scaleFactor), sourceMode(sm), REFluid(REFluid)
{
    SetName("ComptonSource");
    this->id_ne = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Eterm = this->unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    this->id_ntot = this->unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    
    
    //sourceVec = new real_t[kineticGrid->GetNCells()];
    source = new real_t[kineticGrid->GetNCells()];
    
    this->limit = 1000;
    this->wp = gsl_integration_workspace_alloc(limit);
    this->wpOut = gsl_integration_workspace_alloc(limit);
    
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++){
            for(len_t j=0; j<n2[ir]; j++){
                source[offset + j*n1[ir] + i] = EvaluateSource(ir,i,j);
            }
        }
        offset += n1[ir]*n2[ir];
    }
}

/*bool ComptonSource::GridRebuilt(){
    delete [] sourceVec;
    sourceVec = new real_t[this->grid->GetNCells()];
    
    delete [] source;
    source = new real_t[this->grid->GetNCells()];
    
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++){
            for(len_t j=0; j<n2[ir]; j++){
                source[offset + j*n1[ir] + i] = EvaluateSource(ir,i,j);
            }
        }
        offset += n1[ir]*n2[ir];
    }
    
    return true;    
}*/

/*void ComptonSource::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) {
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++){
            for(len_t j=0; j<n2[ir]; j++){
                sourceVec[offset + j*n1[ir] + i] = GetSourceFunction(ir,i,j);
            }
        }
        offset += n1[ir]*n2[ir];
    }
}*/

real_t ComptonSource::innerIntegrand(real_t Eg, void * params){
    real_t r_e = Constants::r0;
    real_t mc2 = Constants::mc2inEV;
    real_t integratedComptonSpectrum = 5.8844; // Integral of the photon flux spectrum over all Eg (in units of mc2).
    
    real_t z = (log(mc2 * Eg/1e6) + 1.2)/0.8;
    real_t GgEg = 1. / integratedComptonSpectrum * exp(- exp(-z) - z + 1.); // * photonFlux (in EvaluateSoure)
    
    real_t p = *(real_t *) params;
    real_t g = sqrt(p*p + 1.);
    real_t cost = 1. - 1. / Eg * (g - 1.) / (Eg + 1. - g);
    real_t dsdO = r_e*r_e / 2.0 * 1.0 / (Eg*Eg) * (Eg / (Eg + 1. - g) + (Eg + 1. - g) / Eg - 1 + cost*cost); // * (Eg - g)*(Eg - g) 
    //real_t dsdp = 1. / ((Eg - g)*(Eg - g)); 
    
    return GgEg * dsdO * mc2; // mc2 from E/mc^2 --> E
}

real_t ComptonSource::integratedPhotonEnergySpectrum(real_t p, void * params){
    if(p == std::numeric_limits<real_t>::infinity())
        return 0;
    real_t g = sqrt(p*p + 1.);
    real_t E0 = (g-1)/2 + sqrt((g-1)*(g-1)/4+(g-1)/2);
    //real_t Emax = inf;
    struct intparams* iparams = reinterpret_cast<struct intparams*>(params);
    
    real_t integral;
    real_t abserr;
    gsl_function F;
    F.function = &(ComptonSource::innerIntegrand);
    F.params = &p;
    gsl_integration_qagiu(&F, E0, 0., 1e-8, iparams->limit, iparams->wp, &integral, &abserr);
    
    real_t S = p / g * integral;
    return S;
}

real_t ComptonSource::integrand(real_t p, void * params){
    if(p == std::numeric_limits<real_t>::infinity())
        return 0;
    real_t S = integratedPhotonEnergySpectrum(p, params);
    real_t p2S = (p*p) * S;
    
    //printf("\np=%.5e, int=%.5e", p, p2S);
    return p2S;
}

/**
 * Evaluates the constant (only grid dependent) source-shape function S(r,p)
 */
real_t ComptonSource::EvaluateSource(len_t ir, len_t i, len_t j) {
    struct intparams params = {this->limit, this->wp};
    if(sourceMode == SOURCE_MODE_FLUID){
        return 0.;
    }
    real_t pm = grid->GetMomentumGrid(ir)->GetP1_f(i);
    real_t pp = grid->GetMomentumGrid(ir)->GetP1_f(i+1);
    if(pm == pp)
        return 0;
    real_t dp = pp-pm;
    real_t pi = (pp+pm)/2.0;
    real_t xim = grid->GetMomentumGrid(ir)->GetP2_f(j);
    real_t xip = grid->GetMomentumGrid(ir)->GetP2_f(j+1);
    real_t dxi = xip-xim;

    real_t integral;
    real_t abserr;
    len_t neval;
    gsl_function F;
    F.function = &(ComptonSource::integrand);
    F.params = &params;
    
    gsl_integration_qng(&F, pm, pp, 0, 1e-8, &integral, &abserr, &neval);
    return scaleFactor * M_PI * this->photonFlux / (dp * dxi * pi*pi) * integral;
}

/**
 * Returns the source at grid point (ir,i,j).
 */
real_t ComptonSource::GetSourceFunction(len_t ir, len_t i, len_t j){
    struct intparams params = {this->limit, this->wp};
    struct intparams paramsOut = {this->limit, this->wpOut};
    if(sourceMode == SOURCE_MODE_FLUID){
        real_t pc = REFluid->GetEffectiveCriticalRunawayMomentum(ir);
        if (isinf(pc))
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
    if(derivId==id_ne) 
        return GetSourceFunction(ir,i,j);
    else if(derivId==id_Eterm) {
        struct intparams params = {this->limit, this->wp};
        real_t pc = REFluid->GetEffectiveCriticalRunawayMomentum(ir);
        if (isinf(pc))
            return 0.;
        real_t Eeff = REFluid->GetEffectiveCriticalField(ir);
        real_t Eterm = unknowns->GetUnknownData(id_Eterm)[ir];
        real_t sgnE = (Eterm>0) - (Eterm<0);
        return -1/2 * (-integrand(pc, &params)) * pc * sgnE/( fabs(Eterm) - Eeff ) ;
    }
    else if(derivId==id_ntot) {
        struct intparams params = {this->limit, this->wp};
        real_t pc = REFluid->GetEffectiveCriticalRunawayMomentum(ir);
        if (isinf(pc))
            return 0.;
        real_t ntot = unknowns->GetUnknownData(id_ntot)[ir];
        return 1/2 * (-integrand(pc, &params)) * pc / ntot;
    }
    else
        return 0;
}

/**
 * Set matrix elements.
 */
/*void ComptonSource::SetMatrixElements(FVM::Matrix *mat, real_t*){
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++){
            for(len_t j=0; j<n2[ir]; j++){
                len_t ind = offset + n1[ir]*j + i;
                mat->SetElement(ind, ir, sourceVec[ind]);
            }
        }
        offset += n1[ir]*n2[ir];
    }
}*/

/**
 * Set vector elements.
 */
/*void ComptonSource::SetVectorElements(real_t *vec, const real_t *x){
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++){
            for(len_t j=0; j<n2[ir]; j++){
                len_t ind = offset + n1[ir]*j + i;
                vec[ind] += sourceVec[ind]*x[ir];
                printf("\nCompton S=%.5e", sourceVec[ind]);
            }
        }
        offset += n1[ir]*n2[ir];
    }
}*/

/**
 * Set jacobian matrix elements.
 */
/*bool ComptonSource::SetJacobianBlock(const len_t , const len_t derivId, FVM::Matrix *jac, const real_t*){
    if(derivId != id_ne)
        return false;
    SetMatrixElements(jac, nullptr);
    return true;
}*/

real_t ComptonSource::EvaluateTotalComptonNumber(real_t pLower, intparams * params, intparams * paramsOut, real_t pUpper){
    real_t integral;
    real_t abserr;
    len_t neval;
    gsl_function F;
    F.function = &(integratedPhotonEnergySpectrum);
    F.params = params;
    if(pUpper == std::numeric_limits<real_t>::infinity()) {
        struct intparams *iparamsOut = reinterpret_cast<struct intparams*>(paramsOut);
        gsl_integration_qagiu(&F, pLower, 0., 1e-8, iparamsOut->limit, iparamsOut->wp, &integral, &abserr);
    } else {
        gsl_integration_qng(&F, pLower, pUpper, 0, 1e-8, &integral, &abserr, &neval);
    }
    return 2. * M_PI * integral;
}
