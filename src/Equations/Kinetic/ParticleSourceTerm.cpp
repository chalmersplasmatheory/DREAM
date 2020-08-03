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
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
        
    // non-trivial temperature jacobian for Maxwellian-shaped particle source
    if(particleSourceShape == PARTICLE_SOURCE_SHAPE_MAXWELLIAN)
        AddUnknownForJacobian(id_Tcold);
}

/**
 * Set the elements in the source function vector.
 */
real_t ParticleSourceTerm::GetSourceFunction(len_t ir, len_t i, len_t j){
    real_t S;
    switch(particleSourceShape){
        case PARTICLE_SOURCE_SHAPE_MAXWELLIAN:{
            real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
            real_t p = grid->GetMomentumGrid(ir)->GetP(i,j);
            S = Constants::RelativisticMaxwellian(p,1,T_cold[ir]);
            break;
        }
        case PARTICLE_SOURCE_SHAPE_DELTA: {
            // XXX: assumes p-xi grid 
            if(i==0)
                S = 1;
            else
                S = 0;
            break;
        }
        default:
            throw FVM::FVMException("ParticleSourceTerm: Invalid particle source shape provided.");
    }
    return S;
}

/**
 * Returns the source function at (ir,i,j) differentiated with respect to the unknown x_derivId at (ir,i,j)
 */
real_t ParticleSourceTerm::GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId){
    real_t dS = 0;
    switch(particleSourceShape){
        case PARTICLE_SOURCE_SHAPE_MAXWELLIAN:
            if(derivId==id_Tcold){
                real_t p = grid->GetMomentumGrid(ir)->GetP(i,j);

                real_t T = unknowns->GetUnknownData(id_Tcold)[ir];
                real_t eps = std::numeric_limits<real_t>::epsilon();
                real_t h = T*sqrt(eps);
                // evaluate numerical temperature derivative
                dS = ( Constants::RelativisticMaxwellian(p,1,T+h ) 
                     - Constants::RelativisticMaxwellian(p,1,T) )
                     / h;
            }
            break;
        case PARTICLE_SOURCE_SHAPE_DELTA: 
            // XXX: assumes p-xi grid 
            dS = 0;
            break;
        default:
            throw FVM::FVMException("ParticleSourceTerm: Invalid particle source shape provided.");
    }
    return dS;
}


