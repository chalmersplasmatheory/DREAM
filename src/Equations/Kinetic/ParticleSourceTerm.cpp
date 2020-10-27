/**
 * Implementation of the particle source term, given by
 *     T = S(r,p) * kappa(r,t),
 * where kappa is the fluid particle source magnitude 
 * which is treated as an unknown quantity, and 
 * S describes the shape of the particle source function
 * in 3-dimensional phase space.
 */

#include "DREAM/Equations/Kinetic/ParticleSourceTerm.hpp"
#include "DREAM/Constants.hpp"
#include <limits>

using namespace DREAM;

/**
 * Constructor
 */
ParticleSourceTerm::ParticleSourceTerm(
    FVM::Grid *kineticGrid, FVM::UnknownQuantityHandler *u, ParticleSourceShape pss
) : FluidSourceTerm(kineticGrid, u), particleSourceShape(pss)
{
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);

    // non-trivial temperature jacobian for Maxwellian-shaped particle source
    if(particleSourceShape == PARTICLE_SOURCE_SHAPE_MAXWELLIAN)
        AddUnknownForJacobian(id_Tcold);
}


/**
 * Normalize the particle source so that it integrates to negative unity
 */
void ParticleSourceTerm::Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *u){
    this->FluidSourceTerm::Rebuild(t,dt,u);
    NormalizeSourceToConstant(-1.0);
}


/**
 * Set the elements in the source function vector.
 */
real_t ParticleSourceTerm::GetSourceFunction(len_t ir, len_t i, len_t j){
    if(grid->IsNegativePitchTrappedIgnorableCell(ir,j))
        return 0;
    switch(particleSourceShape){
        case PARTICLE_SOURCE_SHAPE_MAXWELLIAN:{
            real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
            real_t p = grid->GetMomentumGrid(ir)->GetP(i,j);
            return Constants::RelativisticMaxwellian(p,1,T_cold[ir]);
        }
        case PARTICLE_SOURCE_SHAPE_DELTA:{
            // XXX: assumes p-xi grid 
            len_t n = 2; // only the n innermost p grid points contribute
            if(i<n)
                return n-i;
            else 
                return 0;
        }
        default:
            throw FVM::FVMException("ParticleSourceTerm: Invalid particle source shape provided.");
            return -1;
    }
}

/**
 * Returns the source function at (ir,i,j) differentiated with respect to the unknown x_derivId at (ir,i,j)
 */
real_t ParticleSourceTerm::GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId){
    real_t dS = 0;
    if(grid->IsNegativePitchTrappedIgnorableCell(ir,j))
        return dS;
    switch(particleSourceShape){
        case PARTICLE_SOURCE_SHAPE_MAXWELLIAN:
            if(derivId==id_Tcold){
                real_t p = grid->GetMomentumGrid(ir)->GetP(i,j);
                real_t T = unknowns->GetUnknownData(id_Tcold)[ir];
                Constants::RelativisticMaxwellian(p,1,T,nullptr,&dS);
            }            
            break;
        case PARTICLE_SOURCE_SHAPE_DELTA: 
            break;
        default:
            throw FVM::FVMException("ParticleSourceTerm: Invalid particle source shape provided.");
    }
    return dS;
}


