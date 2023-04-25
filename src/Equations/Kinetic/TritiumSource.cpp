/**
 * Implementation of the kinetic Tritium decay
 * source term, which takes the quadratic form
 *     T = S(p) * n_T(r,t)
 * where S is only a function of phase-space coordinates.
 */

#include "DREAM/Equations/Kinetic/TritiumSource.hpp"
#include "DREAM/Constants.hpp"
//#include <gsl/gsl_integration.h> // ??

using namespace DREAM;

/**
 * Constructor.
 */
TritiumSource::TritiumSource(
    FVM::Grid *kineticGrid, FVM::UnknownQuantityHandler *u, IonHandler *ions, real_t scaleFactor
) : FluidSourceTerm(kineticGrid, u), scaleFactor(scaleFactor)
{
    SetName("TritiumSource");
    this->ions  = ions;
    this->id_nT = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
       
    len_t nT = ions->GetNTritiumIndices();
    sourceVec = new real_t[2 * nT * kineticGrid->GetNCells()];
    
    source = new real_t[kineticGrid->GetNCells()];
    
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++)
            for(len_t j=0; j<n2[ir]; j++){
                source[offset + j*n1[ir] + i] = EvaluateSource(ir,i,j);
            }
        offset += n1[ir]*n2[ir];
    }
}

bool TritiumSource::GridRebuilt(){
    delete [] sourceVec;
    len_t nT = ions->GetNTritiumIndices();
    sourceVec = new real_t[2 * nT * this->grid->GetNCells()];
    
    delete [] source;
    source = new real_t[this->grid->GetNCells()];
    
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++)
            for(len_t j=0; j<n2[ir]; j++){
                source[offset + j*n1[ir] + i] = EvaluateSource(ir,i,j);
            }
        offset += n1[ir]*n2[ir];
    }
    
    return true;    
}

void TritiumSource::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) {
    len_t nT = ions->GetNTritiumIndices();
    const len_t *indT = ions->GetTritiumIndices();
    
    for(len_t iT = 0; iT < nT; iT++){
        len_t offset = indT[iT] * nr * n1[nr-1] * n2[nr-1];
        for(len_t iZ0=0; iZ0<2; iZ0++){
            for(len_t ir=0; ir<nr; ir++){
                for(len_t i=0; i<n1[ir]; i++)
                    for(len_t j=0; j<n2[ir]; j++){
                        sourceVec[offset + j*n1[ir] + i] = GetSourceFunction(ir,i,j);
                    }
                offset += n1[ir]*n2[ir];
            }
        }        
    }
}

real_t TritiumSource::integrand(real_t p, void *){
    real_t Tmax = 18.6e3;
    real_t mc2 = Constants::mc2inEV;
    real_t alpha = 1.0/137.0;
    real_t g = sqrt(p*p + 1);
    real_t fbeta = p * g * (Tmax - mc2 * (g - 1)) * (Tmax - mc2 * (g - 1)) / (1 - exp(-4 * M_PI * alpha * g / p));
    real_t p2fbeta = (p*p) * fbeta;
    return p2fbeta;
}

/**
 * Evaluates the constant (only grid dependent) source-shape function S(r,p)
 */
real_t TritiumSource::EvaluateSource(len_t ir, len_t i, len_t j) {
    real_t C = 1.218e-7; // Normalization factor, 1.2176392e-7
    real_t tau_T = 4500*24*60*60;

    real_t pm = grid->GetMomentumGrid(ir)->GetP1_f(i);
    real_t pp = grid->GetMomentumGrid(ir)->GetP1_f(i+1);
    real_t dp = pp-pm;
    real_t pi = (pp+pm)/2.0;
    
    real_t xim = grid->GetMomentumGrid(ir)->GetP2_f(j);
    real_t xip = grid->GetMomentumGrid(ir)->GetP2_f(j+1);
    real_t dxi = xip-xim;
    
    real_t integral;
    gsl_function F;
    F.function = &(TritiumSource::integrand);
    gsl_integration_qng(&F, pm, pp, 0, 1e-8, &integral, nullptr, nullptr);
    
    return scaleFactor*log(2) / (2.0 * tau_T) * C / (dp * dxi * pi*pi) * integral; 
}

/**
 * Returns the source at grid point (ir,i,j).
 */
real_t TritiumSource::GetSourceFunction(len_t ir, len_t i, len_t j){
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
real_t TritiumSource::GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId){
    if(derivId==id_nT) 
        return GetSourceFunction(ir,i,j);
    else
        return 0;
}

/**
 * Set matrix elements.
 */
void TritiumSource::SetMatrixElements(FVM::Matrix *mat, real_t* /*rhs*/){
    len_t nT = ions->GetNTritiumIndices();
    const len_t *indT = ions->GetTritiumIndices();
    
    for(len_t iT = 0; iT < nT; iT++){
        len_t offset = indT[iT] * nr * n1[nr-1] * n2[nr-1];
        for(len_t iZ0=0; iZ0<2; iZ0++){
            for(len_t ir=0; ir<nr; ir++){
                for(len_t i=0; i<n1[ir]; i++)
                    for(len_t j=0; j<n2[ir]; j++){
                        len_t ind = offset + n1[ir]*j + i;
                        mat->SetElement(ind, ir, sourceVec[ind]);
                    }
                offset += n1[ir]*n2[ir];
            }        
        }
    }
}

/**
 * Set vector elements.
 */
void TritiumSource::SetVectorElements(real_t *vec, const real_t *x){
    len_t nT = ions->GetNTritiumIndices();
    const len_t *indT = ions->GetTritiumIndices();
    
    for(len_t iT = 0; iT < nT; iT++){
        len_t offset = indT[iT] * nr * n1[nr-1] * n2[nr-1];
        for(len_t iZ0=0; iZ0<2; iZ0++){
            for(len_t ir=0; ir<nr; ir++){
                for(len_t i=0; i<n1[ir]; i++)
                    for(len_t j=0; j<n2[ir]; j++){
                        len_t ind = offset + n1[ir]*j + i;
                        vec[ind] += sourceVec[ind]*x[ir];
                    }
                offset += n1[ir]*n2[ir];
            }
        }
    }
}

/**
 * Set jacobian matrix elements.
 */
bool TritiumSource::SetJacobianBlock(const len_t , const len_t derivId, FVM::Matrix *jac, const real_t*){
    if(derivId != id_nT)
        return false;
    SetMatrixElements(jac, nullptr);
    return true;
}
