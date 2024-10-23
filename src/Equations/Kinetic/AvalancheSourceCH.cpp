/**
 * Implementation of the Chiu-Harvey avalanche
 * source term, which takes the quadratic form
 *     T = S(r,p) * n_tot(r,t) * f_re(r,t),
 * where S is only a function of phase-space coordinates.
 */

#include "DREAM/Equations/Kinetic/AvalancheSourceCH.hpp"
#include "DREAM/Constants.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
AvalancheSourceCH::AvalancheSourceCH(
    FVM::Grid *kineticGrid, FVM::UnknownQuantityHandler *u,
    real_t pCutoff, real_t scaleFactor, CHSourcePitchMode sxm
) : FluidSourceTerm(kineticGrid, u), scaleFactor(scaleFactor), sourceXiMode(sxm)
{
    SetName("AvalancheSourceCH");

    id_ntot = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    id_Efield = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    this->pCutoff = pCutoff;
    // non-trivial temperature jacobian for Maxwellian-shaped particle source
    AddUnknownForJacobian(u, id_ntot);
    real_t e = Constants::ec;
    real_t epsmc = Constants::eps0 * Constants::me * Constants::c;
    this->preFactor = (e*e*e*e)/(2*M_PI*epsmc*epsmc*Constants::c);
}

/**
 * Evaluates the constant (only grid dependent) source-shape function S(r,p)
 */
real_t AvalancheSourceCH::EvaluatCHSource(len_t ir, len_t i, len_t j){
    real_t pm = grid->GetMomentumGrid(ir)->GetP1_f(i);
    real_t pp = grid->GetMomentumGrid(ir)->GetP1_f(i+1);
    
    // if pCutoff lies above this cell, return 0.
    // if pCutoff lies inside this cell, use only the corresping fraction 
    // (pp - pCutoff)/(pp - pm) of the source.
    real_t cutFactor = 1.;
    if(pp<=pCutoff)
        return 0;
    else if(pm<pCutoff)
        cutFactor = (pp - pCutoff)/(pp - pm); 
         
    int_t RESign;
    if (this->sourceXiMode == CH_SOURCE_PITCH_ADAPTIVE) {
        const real_t E = unknowns->GetUnknownData(id_Efield)[ir];
        RESign = (E>=0) ? 1: -1;
    } else if (this->sourceXiMode == CH_SOURCE_PITCH_POSITIVE)
        RESign = 1;
    else
        RESign = -1;

    const real_t BA = grid->GetAvalancheCHBounceAverage(ir,i,j, RESign);
    return scaleFactor * preFactor * cutFactor * BA;
}

/**
 * Returns the source at grid point (ir,i,j).
 */
real_t AvalancheSourceCH::GetSourceFunction(len_t ir, len_t i, len_t j){
    real_t S = EvaluateCHSource(ir,i,j);
    const real_t ntot = unknowns->GetUnknownData(id_ntot)[ir];
    return S * ntot;
}

/**
 * Returns the source function at (ir,i,j) differentiated with respect to the unknown x_derivId at (ir,i,j)
 */
real_t AvalancheSourceCH::GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId){
    if(derivId==id_ntot)
        return EvaluateCHSource(ir,i,j);
    else
        return 0;
}


/** 
 * Set matrix elements. 
 */
void AvalancheSourceCH::SetMatrixElements(FVM::Matrix *mat, real_t* /*rhs*/){
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++)
            for(len_t j=0; j<n2[ir]; j++){
                len_t ind = offset + n1[ir]*j + i;
                mat->SetElement(ind, ir, sourceVec[ind] * this->grid->GetMomentumGrid(i)->GetDp2(j) / 2.);
            }
        offset += n1[ir]*n2[ir];
    }        
}


/**
 * Set vector elements.
 */
void AvalancheSourceCH::SetVectorElements(real_t *vec, const real_t *x){
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++)
            for(len_t j=0; j<n2[ir]; j++){
                len_t ind = offset + n1[ir]*j + i;
                real_t f_re_PA = 0;
                for(len_t j_int=0; j_int<n2[ir]; j_int++){
                    f_re_PA += x[offset + n1[ir]*j_int + i] * this->grid->GetMomentumGrid(i)->GetDp2(j_int) / 2.;
                }
                vec[ind] += sourceVec[ind]*f_re_PA;
            }
        offset += n1[ir]*n2[ir];
    }
}


/**
 * Set jacobian matrix elements.
 */
 /* Not needed: the same as in FluidSource
bool AvalancheSourceCH::SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t*){
    if(derivId != uqtyId){
        len_t offset = 0;
        for(len_t ir=0; ir<nr; ir++){
            for(len_t i=0; i<n1[ir]; i++)
                for(len_t j=0; j<n2[ir]; j++){
                    real_t dS = GetSourceFunctionJacobian(ir,i,j,derivId);
                    jac->SetElement(offset + n1[ir]*j + i, ir, dS*x[ir]);
                }
            offset += n1[ir]*n2[ir];
        }
    }
    SetMatrixElements(jac, nullptr);
    return true;
}
*/