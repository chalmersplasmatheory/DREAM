/**
 * Implementation of BounceSurfaceMetric class which is a BounceSurfaceQuantity
 * specified to the metric. The metric is given by the spatial Jacobian (which
 * is described by a FluxSurfaceQuantity) multiplied by the momentum-space jacobian,
 * normalized to a factor of p^2
 */

#include "FVM/Grid/BounceSurfaceMetric.hpp"


using namespace DREAM::FVM;

/** 
 * Constructor
 */ 
BounceSurfaceMetric::BounceSurfaceMetric(Grid *g, FluxSurfaceQuantity *Jacobian, FluxSurfaceQuantity *BOverBmin, FluxSurfaceAverager *FSA)
    : BounceSurfaceQuantity(g,Jacobian), BOverBmin(BOverBmin), fluxSurfaceAverager(FSA){}

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

    bool isPXiGrid = true;
    if(isPXiGrid) {
        FVM::MomentumGrid *mg = grid->GetMomentumGrid(0);
        const real_t *xi0;
        if(fluxGridType == FLUXGRIDTYPE_P2)
            xi0 = mg->GetP2_f();
        else 
            xi0 = mg->GetP2();

        real_t *tmp = new real_t[ntheta_interp_trapped];
        for(len_t ir = 0; ir<nr; ir++)
            for(len_t j=0; j<n2; j++)
                if(IsTrapped(ir,0,j,fluxGridType,grid)){
                    for(len_t it=0; it<ntheta_interp_trapped; it++){
                        real_t theta = ThetaBounceAtIt(ir,0,j,it,fluxGridType);
                        real_t B,Jacobian,ROverR0,NablaR2;
                        fluxSurfaceAverager->GeometricQuantitiesAtTheta(ir,theta,B,Jacobian,ROverR0,NablaR2,fluxGridType);
                        real_t BOverBmin = 1.0;
                        if(Bmin[ir]!=0)
                            BOverBmin = B/Bmin[ir];
                        real_t xiOverXi0 = MomentumGrid::evaluateXiOverXi0(xi0[j],BOverBmin);
                        tmp[it] = Jacobian * mg->evaluatePXiMetricOverP2(xiOverXi0,BOverBmin);
                    }
                    for(len_t i=0; i<n1; i++){ 
                        len_t pind = n1*j+i;
                        bounceData[ir][pind] = new real_t[ntheta_interp_trapped];
                        for(len_t it=0; it<ntheta_interp_trapped; it++)
                            bounceData[ir][pind][it] = tmp[it]; 
                    }
                }
        delete [] tmp;
    } else {
        for(len_t ir = 0; ir<nr; ir++)
            for(len_t i=0; i<n1; i++)
                for(len_t j=0; j<n2; j++)
                    if(IsTrapped(ir,i,j,fluxGridType,grid)){
                        len_t pind = n1*j+i;
                        bounceData[ir][pind] = new real_t[ntheta_interp_trapped];
                        for(len_t it=0; it<ntheta_interp_trapped; it++){
                            real_t theta = ThetaBounceAtIt(ir,i,j,it,fluxGridType);
                            real_t B,Jacobian,ROverR0,NablaR2;
                            fluxSurfaceAverager->GeometricQuantitiesAtTheta(ir,theta,B,Jacobian,ROverR0,NablaR2,fluxGridType);
                            real_t BOverBmin = 1;
                            if(Bmin[ir]!=0)
                                BOverBmin = B/Bmin[ir];

                            grid->GetMomentumGrid(0)->EvaluateMetricOverP2(i,j,fluxGridType, 1, &theta,&BOverBmin, sqrtg_tmp);
                            bounceData[ir][pind][it] = sqrtg_tmp[0] * Jacobian;
                        }
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

    for(len_t ir = 0; ir<nr; ir++){
        const real_t *BOverBmin = this->BOverBmin->GetData(ir,fluxGridType);
        const real_t *Jacobian  = this->fluxSurfaceQuantity->GetData(ir,fluxGridType);
        for(len_t i=0; i<n1; i++)
            for(len_t j=0; j<n2; j++)
                if(!IsTrapped(ir,i,j,fluxGridType,grid)){
                    len_t pind = n1*j+i;
                    bounceData[ir][pind] = new real_t[ntheta_interp_passing];
                    grid->GetMomentumGrid(0)->EvaluateMetricOverP2(i,j,fluxGridType, ntheta_interp_passing, theta_passing, BOverBmin, bounceData[ir][pind]);    
                    for(len_t it=0; it<ntheta_interp_passing; it++)
                        bounceData[ir][pind][it] *= Jacobian[it];
                }
    }
}

/**
 * Evaluates the metric at poloidal angle theta: Jacobian(ir,theta) * sqrtg( B(ir,theta),i,j)
 */
const real_t BounceSurfaceMetric::evaluateAtTheta(len_t ir, len_t i, len_t j, real_t theta, fluxGridType fluxGridType) const {
    real_t ct = cos(theta);
    real_t st = sin(theta);

    return evaluateAtTheta(ir,i,j,theta,ct,st,fluxGridType);
}
/**
 * Evaluates the metric at poloidal angle theta: Jacobian(ir,theta) * sqrtg( B(ir,theta),i,j)
 */
const real_t BounceSurfaceMetric::evaluateAtTheta(len_t /*ir*/, len_t i, len_t j, real_t theta, real_t BOverBmin, real_t Jacobian, fluxGridType fluxGridType) const {
    real_t *sqrtg;
    grid->GetMomentumGrid(0)->EvaluateMetricOverP2(i,j,fluxGridType, 1, &theta,&BOverBmin, sqrtg);
    return *(sqrtg) * Jacobian;
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
    
