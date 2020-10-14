/**
 * Implementation of BounceSurfaceMetric class which is a BounceSurfaceQuantity
 * specified to the metric. The metric is given by the spatial Jacobian (which
 * is described by a FluxSurfaceQuantity) multiplied by the momentum-space jacobian.
 */

#include "FVM/Grid/BounceSurfaceMetric.hpp"


using namespace DREAM::FVM;

/** 
 * Constructor
 */ 
BounceSurfaceMetric::BounceSurfaceMetric(Grid *g, FluxSurfaceQuantity *Jacobian, FluxSurfaceQuantity *B, FluxSurfaceAverager *FSA)
    : BounceSurfaceQuantity(g,Jacobian), B(B), fluxSurfaceAverager(FSA){}

/**
 * Destructor
 */
BounceSurfaceMetric::~BounceSurfaceMetric(){}

/**
 * Set data on passing grid (unlike BounceSurfaceQuantities where this is done 
 * in the FluxSurfaceQuantity instead).
 */
void BounceSurfaceMetric::SetDataForPassing(
    len_t ntheta_interp_passing, const real_t *theta_passing)
{
    InterpolateToFluxGrid(bounceData, FLUXGRIDTYPE_DISTRIBUTION,
        ntheta_interp_passing, theta_passing);
    InterpolateToFluxGrid(bounceData_fr, FLUXGRIDTYPE_RADIAL,
        ntheta_interp_passing, theta_passing);
    InterpolateToFluxGrid(bounceData_f1, FLUXGRIDTYPE_P1,
        ntheta_interp_passing, theta_passing);
    InterpolateToFluxGrid(bounceData_f2, FLUXGRIDTYPE_P2,
        ntheta_interp_passing, theta_passing);

    passingAllocated = true;
}




/**
 * Interpolates data to the poloidal bouncegrid (trapped theta grid).
 * XXX: Assumes same momentum grid at all radii
 */
void BounceSurfaceMetric::InterpolateToBounceGrid(
    real_t ***&bounceData, fluxGridType fluxGridType)
{
    len_t nr = this->nr + (fluxGridType==FLUXGRIDTYPE_RADIAL);
    len_t n1 = np1[0] + (fluxGridType==FLUXGRIDTYPE_P1);
    len_t n2 = np2[0] + (fluxGridType==FLUXGRIDTYPE_P2);
    const real_t *Bmin;
    if(fluxGridType == FLUXGRIDTYPE_RADIAL)
        Bmin = grid->GetRadialGrid()->GetBmin_f();
    else
        Bmin = grid->GetRadialGrid()->GetBmin();

    real_t *sqrtg_tmp = new real_t[1];

    for(len_t ir = 0; ir<nr; ir++)
        for(len_t i=0; i<n1; i++)
            for(len_t j=0; j<n2; j++)
                if(IsTrapped(ir,i,j,fluxGridType,grid)){
                    bounceData[ir][n1*j+i] = new real_t[ntheta_interp_trapped];
                    for(len_t it=0; it<ntheta_interp_trapped; it++){
                        real_t theta = ThetaBounceAtIt(ir,i,j,it,fluxGridType);
                        real_t ct = cos(theta);
                        real_t st = sin(theta);
                        real_t B = fluxSurfaceAverager->BAtTheta(ir, theta,ct,st,fluxGridType);
                        real_t J = fluxSurfaceAverager->JacobianAtTheta(ir, theta,ct,st,fluxGridType);
                        grid->GetMomentumGrid(0)->EvaluateMetric(i,j,fluxGridType, 1, &theta,&B,Bmin[ir], sqrtg_tmp);
                        bounceData[ir][n1*j+i][it] = sqrtg_tmp[0] * J;
                    }
                }
    delete [] sqrtg_tmp;
    trappedAllocated = true;
}

/**
 * Returns stored data.
 */
const real_t *BounceSurfaceMetric::GetData(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const { 
    return GetBounceData(ir,i,j,fluxGridType); 
}

/**
 * Interpolates data to the poloidal fluxgrid theta.
 * XXX: Assumes same momentum grid at all radii
 */
void BounceSurfaceMetric::InterpolateToFluxGrid(
    real_t ***&bounceData, fluxGridType fluxGridType,
    len_t ntheta_interp_passing, const real_t *theta_passing){
    
    len_t nr = this->nr + (fluxGridType==FLUXGRIDTYPE_RADIAL);
    len_t n1 = np1[0] + (fluxGridType==FLUXGRIDTYPE_P1);
    len_t n2 = np2[0] + (fluxGridType==FLUXGRIDTYPE_P2);
    const real_t *Bmin;
    if(fluxGridType == FLUXGRIDTYPE_RADIAL)
        Bmin = grid->GetRadialGrid()->GetBmin_f();
    else
        Bmin = grid->GetRadialGrid()->GetBmin();

    for(len_t ir = 0; ir<nr; ir++){
        for(len_t i=0; i<n1; i++)
            for(len_t j=0; j<n2; j++){
                if(!IsTrapped(ir,i,j,fluxGridType,grid)){
                    len_t pind = n1*j+i;
                    bounceData[ir][pind] = new real_t[ntheta_interp_passing];
                    const real_t *B = this->B->GetData(ir,fluxGridType);
                    grid->GetMomentumGrid(0)->EvaluateMetric(i,j,fluxGridType, ntheta_interp_passing, theta_passing,B,Bmin[ir],bounceData[ir][pind]);    
                    const real_t *Jacobian = this->fluxSurfaceQuantity->GetData(ir,fluxGridType);
                    for(len_t it=0; it<ntheta_interp_passing; it++)
                        bounceData[ir][pind][it] *= Jacobian[it];
                }
            }
    }
}

/**
 * Evaluates the metric at poloidal angle theta: Jacobian(ir,theta) * sqrtg( B(ir,theta),i,j)
 */
const real_t BounceSurfaceMetric::evaluateAtTheta(len_t ir, len_t i, len_t j, real_t theta, fluxGridType fluxGridType) const {
    real_t ct = cos(theta);
    real_t st = sin(theta);
    evaluateAtTheta(ir,i,j,theta,ct,st,fluxGridType);
}
/**
 * Evaluates the metric at poloidal angle theta: Jacobian(ir,theta) * sqrtg( B(ir,theta),i,j)
 */
const real_t BounceSurfaceMetric::evaluateAtTheta(len_t ir, len_t i, len_t j, real_t theta, real_t ct, real_t st, fluxGridType fluxGridType) const {
    real_t Bmin;
    if(fluxGridType == FLUXGRIDTYPE_RADIAL)
        Bmin = grid->GetRadialGrid()->GetBmin_f(ir);
    else
        Bmin = grid->GetRadialGrid()->GetBmin(ir);
    real_t B = fluxSurfaceAverager->BAtTheta(ir,theta,ct,st,fluxGridType);
    real_t J = fluxSurfaceAverager->JacobianAtTheta(ir,theta,ct,st,fluxGridType);
    real_t *sqrtg;
    grid->GetMomentumGrid(0)->EvaluateMetric(i,j,fluxGridType, 1, &theta,&B,Bmin, sqrtg);
    return *(sqrtg) * J;
}

/**
 * Deallocator
 */
void BounceSurfaceMetric::DeleteData(real_t ***&data, bool **isTrapped, len_t nr, len_t np1, len_t np2){
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i = 0; i<np1*np2; i++){
            bool trap = isTrapped[ir][i];
            // delete only if allocated
            if( (trap && trappedAllocated) || ( (!trap) && passingAllocated) ) 
                delete [] data[ir][i];   
        }         
        delete [] data[ir];
    }
    delete [] data;
}
    
