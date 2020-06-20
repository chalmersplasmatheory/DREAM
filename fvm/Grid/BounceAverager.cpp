/**
 * Implementation of the FluxSurfaceAverager class which 
 * handles everything needed to carry out flux surface 
 * averages in DREAM. Is initialized by a RadialGridGenerator
 * via the method SetReferenceMagneticFieldData.
 */

#include "FVM/Grid/BounceAverager.hpp"

using namespace std;
using namespace DREAM::FVM;

/**
 * Constructor.
 */
BounceAverager::BounceAverager(
    Grid *g, FluxSurfaceAverager* fsa, len_t ntheta_interp_trapped,
    enum OptionConstants::momentumgrid_type mgtype,
    FluxSurfaceAverager::quadrature_method q_method_trapped
) : grid(g), fluxSurfaceAverager(fsa), ntheta_interp_trapped(ntheta_interp_trapped) {
    nr = grid->GetNr();
    np1 = new len_t[nr];
    np2 = new len_t[nr];
    for(len_t ir = 0; ir<nr; ir++){
        MomentumGrid *mg = grid->GetMomentumGrid(ir);
        np1[ir] = mg->GetNp1();
        np2[ir] = mg->GetNp2();
    }
    gridModePXI = (mgtype == OptionConstants::MOMENTUMGRID_TYPE_PXI);

    geometryIsSymmetric = fluxSurfaceAverager->isGeometrySymmetric();
    ntheta_interp       = fluxSurfaceAverager->GetNTheta();
    theta               = fluxSurfaceAverager->GetTheta();
    weights             = fluxSurfaceAverager->GetWeights();
    integratePassingAdaptive = fluxSurfaceAverager->isIntegrationAdaptive();

    InitializeQuadrature(q_method_trapped);

    gsl_acc = gsl_interp_accel_alloc();
}

BounceAverager::~BounceAverager(){
    gsl_interp_accel_free(gsl_acc);
    if(integrateTrappedAdaptive)
        gsl_integration_workspace_free(gsl_adaptive);
    else
        gsl_integration_fixed_free(gsl_w);
    delete [] np1;
    delete [] np2;

    DeallocateBounceIntegralQuantities();
    DeallocateBounceGridMagneticQuantities(
        theta_bounceGrid, weights_bounceGrid,
        B_bounceGrid, ROverR0_bounceGrid, 
        NablaR2_bounceGrid, FLUXGRIDTYPE_DISTRIBUTION);
    DeallocateBounceGridMagneticQuantities(
        theta_bounceGrid_fr, weights_bounceGrid_fr,
        B_bounceGrid_fr, ROverR0_bounceGrid_fr, 
        NablaR2_bounceGrid_fr, FLUXGRIDTYPE_RADIAL);
    DeallocateBounceGridMagneticQuantities(
        theta_bounceGrid_f1, weights_bounceGrid_f1,
        B_bounceGrid_f1, ROverR0_bounceGrid_f1, 
        NablaR2_bounceGrid_f1, FLUXGRIDTYPE_P1);
    DeallocateBounceGridMagneticQuantities(
        theta_bounceGrid_f2, weights_bounceGrid_f2,
        B_bounceGrid_f2, ROverR0_bounceGrid_f2, 
        NablaR2_bounceGrid_f2, FLUXGRIDTYPE_P2);


}



// Evaluates the bounce average {F} of a function F = F(xi/xi0, B/Bmin, R/R0) on grid point (ir,i,j). 
real_t BounceAverager::CalculateBounceAverage(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F){
    /**
     * XXX: Here we assume same grid at all radii.
     */
    MomentumGrid *mg = grid->GetMomentumGrid(0);
    real_t Bmin = fluxSurfaceAverager->GetBmin(ir,fluxGridType == FLUXGRIDTYPE_RADIAL);
    real_t Bmax = fluxSurfaceAverager->GetBmax(ir,fluxGridType == FLUXGRIDTYPE_RADIAL);

    real_t xi0;
    if (fluxGridType == FLUXGRIDTYPE_P1){
        xi0 = mg->GetXi0_f1(i,j);
    } else if (fluxGridType == FLUXGRIDTYPE_P2){
        xi0 = mg->GetXi0_f2(i,j);
    } else
        xi0 = mg->GetXi0(i,j);
    
    if(xi0*xi0 < 1e-30){
        real_t BounceAverageSingular, sqrtterm;
        for(len_t it=0; it<ntheta_interp; it++){
            sqrtterm = sqrt(1-quad_x_ref[it]*quad_x_ref[it]);
            BounceAverageSingular +=  quad_w_ref[it]* (F(sqrtterm,1,1,1)+F(-sqrtterm,1,1,1))/sqrtterm;
        }
        BounceAverageSingular /= M_PI;
        return BounceAverageSingular;
    } else
        return EvaluateBounceIntegral(ir,i,j,fluxGridType, F) / GetVp(ir, i, j,fluxGridType);
}




const real_t* BounceAverager::GetB(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const {
    len_t n1 = np1[ir]+(fluxGridType==FLUXGRIDTYPE_P1);
    len_t ind = j*n1+i;
        
    if (GetIsTrapped(ir,i,j,fluxGridType)) {
        switch (fluxGridType){
            case FLUXGRIDTYPE_DISTRIBUTION:
                return B_bounceGrid[ir][ind];
            case FLUXGRIDTYPE_RADIAL:
                return B_bounceGrid_fr[ir][ind];
            case FLUXGRIDTYPE_P1:
                return B_bounceGrid_f1[ir][ind+j];
            case FLUXGRIDTYPE_P2:
                return B_bounceGrid_f2[ir][ind];
            default: 
                return nullptr;
        }
    } else {
        return fluxSurfaceAverager->GetB(ir, fluxGridType==FLUXGRIDTYPE_RADIAL);
    }
}

const real_t* BounceAverager::GetROverR0(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const {
    len_t n1 = np1[ir]+(fluxGridType==FLUXGRIDTYPE_P1);
    len_t ind = j*n1+i;
    
    if (GetIsTrapped(ir,i,j,fluxGridType)) {
        switch (fluxGridType){
            case FLUXGRIDTYPE_DISTRIBUTION:
                return ROverR0_bounceGrid[ir][ind];
            case FLUXGRIDTYPE_RADIAL:
                return ROverR0_bounceGrid_fr[ir][ind];
            case FLUXGRIDTYPE_P1:
                return ROverR0_bounceGrid_f1[ir][ind];
            case FLUXGRIDTYPE_P2:
                return ROverR0_bounceGrid_f2[ir][ind];
            default: 
                throw FVMException("Invalid fluxGridType: '%d' called in BounceAverager getter.", fluxGridType);
        }
    } else {
        return fluxSurfaceAverager->GetROverR0(ir, fluxGridType==FLUXGRIDTYPE_RADIAL);
    }
}

const real_t* BounceAverager::GetNablaR2(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const {
    len_t n1 = np1[ir]+(fluxGridType==FLUXGRIDTYPE_P1);
    len_t ind = j*n1+i;
    
    if (GetIsTrapped(ir,i,j,fluxGridType)) {
        switch (fluxGridType){
            case FLUXGRIDTYPE_DISTRIBUTION:
                return NablaR2_bounceGrid[ir][ind];
            case FLUXGRIDTYPE_RADIAL:
                return NablaR2_bounceGrid_fr[ir][ind];
            case FLUXGRIDTYPE_P1:
                return NablaR2_bounceGrid_f1[ir][ind];
            case FLUXGRIDTYPE_P2:
                return NablaR2_bounceGrid_f2[ir][ind];
            default: 
                throw FVMException("Invalid fluxGridType: '%d' called in BounceAverager getter.", fluxGridType);
        }
    } else {
        return fluxSurfaceAverager->GetNablaR2(ir, fluxGridType==FLUXGRIDTYPE_RADIAL);
    }
}

const real_t* BounceAverager::GetTheta(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const {
    len_t n1 = np1[ir]+(fluxGridType==FLUXGRIDTYPE_P1);
    len_t ind = j*n1+i;
    
    if (GetIsTrapped(ir,i,j,fluxGridType)) {
        switch (fluxGridType){
            case FLUXGRIDTYPE_DISTRIBUTION:
                return theta_bounceGrid[ir][ind];
            case FLUXGRIDTYPE_RADIAL:
                return theta_bounceGrid_fr[ir][ind];
            case FLUXGRIDTYPE_P1:
                return theta_bounceGrid_f1[ir][ind+j];
            case FLUXGRIDTYPE_P2:
                return theta_bounceGrid_f2[ir][ind];
            default: 
                throw FVMException("Invalid fluxGridType: '%d' called in BounceAverager getter.", fluxGridType);
        }
    } else {
        return fluxSurfaceAverager->GetTheta();
    }
}
const real_t* BounceAverager::GetWeights(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const {
    len_t n1 = np1[ir]+(fluxGridType==FLUXGRIDTYPE_P1);
    len_t ind = j*n1+i;
    if (GetIsTrapped(ir,i,j,fluxGridType)) {
        switch (fluxGridType){
            case FLUXGRIDTYPE_DISTRIBUTION:
                return weights_bounceGrid[ir][ind];
            case FLUXGRIDTYPE_RADIAL:
                return weights_bounceGrid_fr[ir][ind];
            case FLUXGRIDTYPE_P1:
                return weights_bounceGrid_f1[ir][ind+j];
            case FLUXGRIDTYPE_P2:
                return weights_bounceGrid_f2[ir][ind];
            default: 
                throw FVMException("Invalid fluxGridType: '%d' called in BounceAverager getter.", fluxGridType);
        }
    } else 
        return fluxSurfaceAverager->GetWeights();       
}

const real_t* BounceAverager::GetMetric(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const {
    len_t n1 = np1[ir]+(fluxGridType==FLUXGRIDTYPE_P1);
    len_t ind = j*n1+i;
    switch (fluxGridType){
        case FLUXGRIDTYPE_DISTRIBUTION:
            return metricSqrtG[ir][ind];
        case FLUXGRIDTYPE_RADIAL:
            return metricSqrtG_fr[ir][ind];
        case FLUXGRIDTYPE_P1:
            return metricSqrtG_f1[ir][ind+j];
        case FLUXGRIDTYPE_P2:
            return metricSqrtG_f2[ir][ind];
        default: 
            throw FVMException("Invalid fluxGridType: '%d' called in BounceAverager getter.", fluxGridType);
    }
}
const real_t BounceAverager::GetVp(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const{
    len_t n1 = np1[ir];
    len_t n2 = np2[ir];
    if(fluxGridType == FLUXGRIDTYPE_DISTRIBUTION)
        return Vp[ir][n1*j+i];    
    else if(fluxGridType == FLUXGRIDTYPE_RADIAL)
        return Vp_fr[ir][n1*j+i];
    else if(fluxGridType == FLUXGRIDTYPE_P1)
        return Vp_f1[ir][(n1+1)*j+i];
    else if(fluxGridType == FLUXGRIDTYPE_P2)
        return Vp_f2[ir][n1*j+i];
    else
        throw FVMException("Invalid fluxGridType: '%d' called in BounceAverager getter.", fluxGridType);

} 


const real_t BounceAverager::GetIsTrapped(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const{
    len_t n1 = np1[ir];
    len_t n2 = np2[ir];
    if(fluxGridType == FLUXGRIDTYPE_DISTRIBUTION)
        return isTrapped[ir][n1*j+i];    
    else if(fluxGridType == FLUXGRIDTYPE_RADIAL)
        return isTrapped_fr[ir][n1*j+i];
    else if(fluxGridType == FLUXGRIDTYPE_P1)
        return isTrapped_f1[ir][(n1+1)*j+i];
    else if(fluxGridType == FLUXGRIDTYPE_P2)
        return isTrapped_f2[ir][n1*j+i];
    else
        throw FVMException("Invalid fluxGridType: '%d' called in BounceAverager getter.", fluxGridType);
} 



real_t BounceAverager::EvaluateBounceIntegral(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F){
        /**
     * XXX: Here we assume same grid at all radii.
     */
    MomentumGrid *mg = grid->GetMomentumGrid(0);
    real_t Bmin = fluxSurfaceAverager->GetBmin(ir,fluxGridType == FLUXGRIDTYPE_RADIAL);
    real_t Bmax = fluxSurfaceAverager->GetBmax(ir,fluxGridType == FLUXGRIDTYPE_RADIAL);

    real_t xi0;
    if (fluxGridType == FLUXGRIDTYPE_P1){
        xi0 = mg->GetXi0_f1(i,j);
    } else if (fluxGridType == FLUXGRIDTYPE_P2){
        xi0 = mg->GetXi0_f2(i,j);
    } else
        xi0 = mg->GetXi0(i,j);

    std::function<real_t(real_t,real_t,real_t,real_t)> F_eff;
    // If trapped, adds contribution from -xi0, since negative xi0 are presumably not kept on the grid.
    if (GetIsTrapped(ir,i,j,fluxGridType))
        F_eff = [&](real_t x, real_t  y, real_t z, real_t w){return  F(x,y,z,w) + F(-x,y,z,w) ;};
    else 
        F_eff = F;

    const real_t *B       = GetB(ir,i,j,fluxGridType);
    const real_t *ROverR0 = GetROverR0(ir,i,j,fluxGridType);
    const real_t *NablaR2 = GetNablaR2(ir,i,j,fluxGridType);
    const real_t *weights = GetWeights(ir,i,j,fluxGridType);
    const real_t *sqrtg   = GetMetric(ir,i,j,fluxGridType);

    if(xi0*xi0 < 1e-30)
        return 0;

    real_t xiOverXi0,w,BOverBmin;        
    real_t BounceIntegral = 0;
    for (len_t it = 0; it<ntheta_interp; it++) {
        real_t xi2 = 1- B[it]/Bmin * (1-xi0*xi0);
        xiOverXi0 = sqrt(xi2/(xi0*xi0));
        w = weights[it];
        BOverBmin = B[it]/Bmin;

        BounceIntegral += 2*M_PI*w*sqrtg[it]*F_eff(xiOverXi0,BOverBmin,ROverR0[it],NablaR2[it]);
    }        
    return BounceIntegral;
    
}

// Initializes quadrature quantities. See FluxSurfaceAverager implementation
// for more details.
void BounceAverager::InitializeQuadrature(FluxSurfaceAverager::quadrature_method q_method){
    const gsl_integration_fixed_type *quadratureRule;
    std::function<real_t(real_t,real_t,real_t)> QuadFunc;
    switch(q_method){
        case FluxSurfaceAverager::QUAD_FIXED_LEGENDRE:
            // Legendre quadrature is often good for smooth finite functions on a finite interval
            quadratureRule = gsl_integration_fixed_legendre;
            QuadFunc = [](real_t /*x*/, real_t /*x_min*/, real_t /*x_max*/)
                                    {return 1;};
            //quadratureWeightFunction = [](real_t /*x*/, real_t /*x_min*/, real_t /*x_max*/)
            //                        {return 1;};
            break;
        case FluxSurfaceAverager::QUAD_FIXED_CHEBYSHEV:
            // Chebyshev quadrature may be better for integration along trapped orbits
            // where the metric takes the form of this weight function.
            quadratureRule = gsl_integration_fixed_chebyshev;
            QuadFunc = [](real_t x, real_t x_min, real_t x_max)
                                    {return 1/sqrt((x_max-x)*(x-x_min) );};
            //quadratureWeightFunction = [](real_t x, real_t x_min, real_t x_max)
            //                        {return 1/sqrt((x_max-x)*(x-x_min) );};
            break;
        case FluxSurfaceAverager::QUAD_ADAPTIVE:
            gsl_adaptive = gsl_integration_workspace_alloc(1000);
            integrateTrappedAdaptive = true;
            return;
        default:
            throw FVMException("Quadrature rule '%d' not supported by FluxSurfaceAverager.", q_method);            
    }
    real_t unusedArgument = 0;
    gsl_w = gsl_integration_fixed_alloc(quadratureRule,ntheta_interp_trapped,0,1,unusedArgument,unusedArgument);
    quad_x_ref = gsl_w->x;
    quad_w_ref = gsl_w->weights;
    for(len_t it=0; it<ntheta_interp_trapped; it++)
        quad_w_ref[it]/= QuadFunc(quad_x_ref[it],0,1);

}


void BounceAverager::Rebuild(){
    bool hasTrapped = InitializeBounceIntegralQuantities();

    // If using fixed quadrature on passing grid, store metric on theta grid
    if(!integratePassingAdaptive){
        SetMetricOnPassingGrid(metricSqrtG, isTrapped, FLUXGRIDTYPE_DISTRIBUTION);
        SetMetricOnPassingGrid(metricSqrtG_fr, isTrapped_fr, FLUXGRIDTYPE_RADIAL);
        SetMetricOnPassingGrid(metricSqrtG_f1, isTrapped_f1, FLUXGRIDTYPE_P1);
        SetMetricOnPassingGrid(metricSqrtG_f2, isTrapped_f2, FLUXGRIDTYPE_P2);
    }

    // If no particles are trapped or trapped orbits are 
    // integrated adaptively, this is all we need
    if( (!hasTrapped) ||  integrateTrappedAdaptive)
        return;


    // Otherwise, for fixed quadratures, we interpolate 
    // everything to the quadrature grid
    InterpolateMagneticQuantitiesToBounceGrid(
        theta_bounceGrid, weights_bounceGrid, theta_b1, 
        theta_b2, B_bounceGrid, ROverR0_bounceGrid, 
        NablaR2_bounceGrid, metricSqrtG, FLUXGRIDTYPE_DISTRIBUTION);
    InterpolateMagneticQuantitiesToBounceGrid(
        theta_bounceGrid_fr, weights_bounceGrid_fr, theta_b1_fr,  
        theta_b2_fr, B_bounceGrid_fr, ROverR0_bounceGrid_fr, 
        NablaR2_bounceGrid_fr, metricSqrtG_fr, FLUXGRIDTYPE_RADIAL);
    InterpolateMagneticQuantitiesToBounceGrid(
        theta_bounceGrid_f1, weights_bounceGrid_f1, theta_b1_f1, 
        theta_b2_f1, B_bounceGrid_f1, ROverR0_bounceGrid_f1, 
        NablaR2_bounceGrid_f1, metricSqrtG_f1, FLUXGRIDTYPE_P1);
    InterpolateMagneticQuantitiesToBounceGrid(
        theta_bounceGrid_f2, weights_bounceGrid_f2, theta_b1_f2, 
        theta_b2_f2, B_bounceGrid_f2, ROverR0_bounceGrid_f2, 
        NablaR2_bounceGrid_f2, metricSqrtG_f2, FLUXGRIDTYPE_P2);


    SetVp(Vp,FLUXGRIDTYPE_DISTRIBUTION);
    SetVp(Vp_fr,FLUXGRIDTYPE_RADIAL);
    SetVp(Vp_f1,FLUXGRIDTYPE_P1);
    SetVp(Vp_f2,FLUXGRIDTYPE_P2);

    /**
     * Do something like 
     * grid->SetBounceIntegralQuantities(
     *  Vp,Vp_fr,Vp_f1,Vp_f2,
     *  isTrapped, isTrapped_fr, isTrapped_f1, isTrapped_f2,
     *  theta_b1, theta_b1_fr, theta_b1_f1, theta_b1_f2,
     *  theta_b2, theta_b2_fr, theta_b2_f1, theta_b2_f2
     * );
     */
}

void BounceAverager::SetVp(real_t**&Vp, fluxGridType fluxGridType){
    function<real_t(real_t,real_t,real_t,real_t)> unityFunc 
               = [](real_t,real_t,real_t,real_t){return 1;};

    len_t nr = this->nr + (fluxGridType == FLUXGRIDTYPE_RADIAL);
    len_t n1, n2;
    for(len_t ir = 0; ir<nr; ir++){
        n1 = np1[ir] + (fluxGridType == FLUXGRIDTYPE_P1);
        n2 = np2[ir] + (fluxGridType == FLUXGRIDTYPE_P2);
        for(len_t j = 0; j<n2; j++){
            for(len_t i = 0; i<n1; i++){
                len_t pind = j*n1+i;
                Vp[ir][pind] = EvaluateBounceIntegral(ir,i,j,fluxGridType,unityFunc);
            }
        }
    }

}

void BounceAverager::SetMetricOnPassingGrid(real_t ***&metricSqrtG, bool **isTrapped, fluxGridType fluxGridType){
    real_t Bmin;
    bool rFluxGrid = (fluxGridType == FLUXGRIDTYPE_RADIAL);
    len_t nr = this->nr + rFluxGrid;
    MomentumGrid *mg;
    len_t n1, n2;
    for(len_t ir = 0; ir<nr; ir++){
        /**
         * XXX: Here we assume same grid at all radii.
         */
        mg = grid->GetMomentumGrid(0);
        Bmin = fluxSurfaceAverager->GetBmin(ir,rFluxGrid);

        n1 = np1[ir] + (fluxGridType == FLUXGRIDTYPE_P1);
        n2 = np2[ir] + (fluxGridType == FLUXGRIDTYPE_P2);

        const real_t *B = fluxSurfaceAverager->GetB(ir,rFluxGrid);
        const real_t *Jacobian = fluxSurfaceAverager->GetJacobian(ir,rFluxGrid);

        for(len_t j = 0; j<n2; j++){
            for(len_t i = 0; i<n1; i++){
                len_t pind = j*n1+i;
                if(!isTrapped[ir][pind]){
                    mg->EvaluateMetric(i,j,fluxGridType, ntheta_interp,  
                                theta, B, Bmin, metricSqrtG[ir][pind]);
                    for(len_t it=0; it<ntheta_interp; it++)
                        metricSqrtG[ir][pind][it] *= Jacobian[it];
                }
            }
        }
    }
}

void BounceAverager::InterpolateMagneticQuantitiesToBounceGrid(
    real_t ***&theta_bounceGrid, real_t ***&weights_bounceGrid, real_t **theta_b1, real_t **theta_b2, 
    real_t ***&B_bounceGrid, real_t ***&ROverR0_bounceGrid, 
    real_t ***&NablaR2_bounceGrid, real_t ***&metricSqrtG, fluxGridType fluxGridType
){
    AllocateBounceGridMagneticQuantities(
        theta_bounceGrid, weights_bounceGrid,
        B_bounceGrid, ROverR0_bounceGrid, 
        NablaR2_bounceGrid, fluxGridType);

    real_t Bmin;
    bool rFluxGrid = (fluxGridType == FLUXGRIDTYPE_RADIAL);
    len_t nr = this->nr + rFluxGrid;
    MomentumGrid *mg;
    len_t n1, n2;
    for(len_t ir = 0; ir<nr; ir++){
        /**
         * XXX: Here we assume same grid at all radii.
         */
        mg = grid->GetMomentumGrid(0);
        Bmin = fluxSurfaceAverager->GetBmin(ir,rFluxGrid);

        n1 = np1[ir] + (fluxGridType == FLUXGRIDTYPE_P1);
        n2 = np2[ir] + (fluxGridType == FLUXGRIDTYPE_P2);

        for(len_t j = 0; j<n2; j++){
            for(len_t i = 0; i<n1; i++){
                len_t pind = j*n1+i;
                if(isTrapped[ir][pind]){
                    for (len_t it=0; it<ntheta_interp_trapped; it++) {
                        theta_bounceGrid[ir][pind][it]    = theta_b1[ir][pind] + (theta_b2[ir][pind]-theta_b1[ir][pind]) * quad_x_ref[it] ;
                        weights_bounceGrid[ir][pind][it]  = (theta_b2[ir][pind] - theta_b1[ir][pind]) * quad_w_ref[it];
                        B_bounceGrid[ir][pind][it]        = fluxSurfaceAverager->evaluateBAtTheta(ir,theta_bounceGrid[ir][pind][it],rFluxGrid);
                        ROverR0_bounceGrid[ir][pind][it]  = fluxSurfaceAverager->evaluateROverR0AtTheta(ir,theta_bounceGrid[ir][pind][it],rFluxGrid);
                        NablaR2_bounceGrid[ir][pind][it]  = fluxSurfaceAverager->evaluateNablaR2AtTheta(ir,theta_bounceGrid[ir][pind][it],rFluxGrid);
                    }
                    mg->EvaluateMetric(i,j,fluxGridType, ntheta_interp_trapped,  
                        theta_bounceGrid[ir][pind], B_bounceGrid[ir][pind], Bmin, metricSqrtG[ir][pind]);
                    for (len_t it=0; it<ntheta_interp_trapped; it++) {
                        metricSqrtG[ir][pind][it] *= fluxSurfaceAverager->evaluateJacobianAtTheta(ir,theta_bounceGrid[ir][pind][it],rFluxGrid);
                    }

                }
            }
        }
    }
}

bool BounceAverager::InitializeBounceIntegralQuantities(){
    AllocateBounceIntegralQuantities();
    bool hasTrapped;
    hasTrapped  = SetIsTrapped(isTrapped,    theta_b1,    theta_b2,    FLUXGRIDTYPE_DISTRIBUTION);
    hasTrapped += SetIsTrapped(isTrapped_fr, theta_b1_fr, theta_b2_fr, FLUXGRIDTYPE_RADIAL);
    hasTrapped += SetIsTrapped(isTrapped_f1, theta_b1_f1, theta_b2_f1, FLUXGRIDTYPE_P1);
    hasTrapped += SetIsTrapped(isTrapped_f2, theta_b1_f2, theta_b2_f2, FLUXGRIDTYPE_P2);

    return hasTrapped;
}


bool BounceAverager::SetIsTrapped(bool **&isTrapped, real_t **&theta_b1, real_t **&theta_b2, fluxGridType fluxGridType){
    bool hasTrapped = false;
    real_t Bmin, Bmax;
    len_t nr = this->nr + (fluxGridType == FLUXGRIDTYPE_RADIAL);
    MomentumGrid *mg;
    len_t n1, n2;
    for(len_t ir = 0; ir<nr; ir++){
        /**
         * XXX: Here we assume same grid at all radii.
         */
        mg = grid->GetMomentumGrid(0);
        Bmin = fluxSurfaceAverager->GetBmin(ir,fluxGridType == FLUXGRIDTYPE_RADIAL);
        Bmax = fluxSurfaceAverager->GetBmax(ir,fluxGridType == FLUXGRIDTYPE_RADIAL);

        n1 = np1[ir];
        n2 = np2[ir];
        const real_t* xi0;
        if (fluxGridType == FLUXGRIDTYPE_P1){
            n1++;
            xi0 = mg->GetXi0_f1();
        } else if (fluxGridType == FLUXGRIDTYPE_P2){
            n2++;
            xi0 = mg->GetXi0_f2();
        } else
            xi0 = mg->GetXi0();
        // in cylindrical grid or at r=0, isTrapped is false and we skip to next radius 
        if(Bmin==Bmax){
            for(len_t i = 0; i<mg->GetNCells(); i++)
                isTrapped[ir][i] = false;
            continue;
        }

        for(len_t j = 0; j<n2; j++){
            for(len_t i = 0; i<n1; i++){
                len_t pind = j*n1+i;
                if((1-xi0[pind]*xi0[pind]) > Bmin/Bmax){
                    isTrapped[ir][pind] = true;
                    hasTrapped = true;
                    FindBouncePoints(ir, xi0[pind], fluxGridType==FLUXGRIDTYPE_RADIAL,
                                    &theta_b1[ir][pind],&theta_b2[ir][pind]);
                } else {
                    isTrapped[ir][n1*j+i] = false;

                }
            }
        }
    }
    return hasTrapped;
}


/**
 * The function is used in the evaluation of the effective passing fraction,
 * and represents x / <1-x B/Bmax>
 */
struct xiFuncParams {real_t xi0; len_t ir; bool rFluxGrid; FluxSurfaceAverager *fluxSurfaceAverager;};
real_t BounceAverager::xiParticleFunction(real_t theta, void *p){
    struct xiFuncParams *params = (struct xiFuncParams *) p;
    
    real_t xi0 = params->xi0; 
    len_t ir   = params->ir;
    bool rFluxGrid = params->rFluxGrid;
    FluxSurfaceAverager *fsa = params->fluxSurfaceAverager;
    return 1 - (1-xi0*xi0) *fsa->evaluateBAtTheta(ir,theta,rFluxGrid) / fsa->GetBmin(ir, rFluxGrid) ;
}

// calculates theta_bounce1 and theta_bounce2 with a root finding algorithm
void BounceAverager::FindBouncePoints(len_t ir, real_t xi0, bool rFluxGrid, real_t *theta_b1, real_t *theta_b2){
    // define GSL function xi_particle as function of theta which evaluates B_interpolator(_fr) at theta.

    xiFuncParams xi_params = {xi0,ir,rFluxGrid,fluxSurfaceAverager}; 
    gsl_function gsl_func;
    gsl_func.function = &(xiParticleFunction);
    gsl_func.params = &xi_params;


    // A guess can be estimated assuming that B ~ 1/[R+rcos(theta)] 
    // -- probably often a good guess. Then:
    // real_t cos_theta_b_guess = ((1-xi0*xi0)/Bmin - (1/Bmin + 1/Bmax)/2 ) * 2/( 1/Bmin - 1/Bmax );

    // Look for first root between 0 and pi
    real_t x_lower = 0.0;
    real_t x_upper = M_PI;
    real_t root = 0;
    FindThetaBounceRoots(&x_lower, &x_upper, &root, gsl_func);
    
    // In symmetric field, this root corresponds to theta_b2
    if (geometryIsSymmetric){
        *theta_b2 = x_lower;
        *theta_b1 = -x_lower;
    } else {

        // if xi(theta = x_upper + epsilon) is real, the root 
        // corresponds to the lower bounce point theta_b1 
        if ( xiParticleFunction(x_upper,&xi_params) > 0 ){
            *theta_b1 = x_upper;
            x_lower = M_PI; 
            x_upper = 2*M_PI;
            FindThetaBounceRoots(&x_lower, &x_upper,&root, gsl_func);
            *theta_b2 = x_lower;
        } else {
            *theta_b2 = x_lower;
            x_lower = -M_PI;
            x_upper = 0;
            FindThetaBounceRoots(&x_lower, &x_upper, &root, gsl_func);
            *theta_b1 = x_upper; 
        }
    }
}


void BounceAverager::DeallocateBounceGridMagneticQuantities(real_t ***&theta_bounceGrid, real_t ***&weights_bounceGrid,  real_t ***&B_bounceGrid, 
            real_t ***&ROverR0_bounceGrid, real_t ***&NablaR2_bounceGrid, fluxGridType fluxGridType){
            
    if(theta_bounceGrid == nullptr)
        return;
    len_t nr = this->nr + (fluxGridType == FLUXGRIDTYPE_RADIAL);
    len_t n1, n2, pind;
    for(len_t ir = 0; ir<nr; ir++){
        n1 = np1[ir] + (fluxGridType == FLUXGRIDTYPE_P1);
        n2 = np2[ir] + (fluxGridType == FLUXGRIDTYPE_P2);
        for(len_t j = 0; j<n2; j++){
            for(len_t i = 0; i<n1; i++){
                pind = j*n1+i;
                delete [] theta_bounceGrid[ir][pind];
                delete [] weights_bounceGrid[ir][pind];
                delete [] B_bounceGrid[ir][pind];
                delete [] ROverR0_bounceGrid[ir][pind];
                delete [] NablaR2_bounceGrid[ir][pind];
            }
        }
        delete [] theta_bounceGrid[ir];
        delete [] weights_bounceGrid[ir];
        delete [] B_bounceGrid[ir];
        delete [] ROverR0_bounceGrid[ir];
        delete [] NablaR2_bounceGrid[ir];
    }
    delete [] theta_bounceGrid;
    delete [] weights_bounceGrid;
    delete [] B_bounceGrid;
    delete [] ROverR0_bounceGrid;
    delete [] NablaR2_bounceGrid;    
}

void BounceAverager::AllocateBounceGridMagneticQuantities(real_t ***&theta_bounceGrid, real_t ***&weights_bounceGrid,  real_t ***&B_bounceGrid, 
            real_t ***&ROverR0_bounceGrid, real_t ***&NablaR2_bounceGrid, fluxGridType fluxGridType){

    DeallocateBounceGridMagneticQuantities(
        theta_bounceGrid, weights_bounceGrid,
        B_bounceGrid, ROverR0_bounceGrid, 
        NablaR2_bounceGrid, fluxGridType);
    
    len_t nr = this->nr + (fluxGridType == FLUXGRIDTYPE_RADIAL);
    len_t n1, n2, pind;

    theta_bounceGrid   = new real_t**[nr];
    weights_bounceGrid = new real_t**[nr];
    B_bounceGrid       = new real_t**[nr];
    ROverR0_bounceGrid = new real_t**[nr];
    NablaR2_bounceGrid = new real_t**[nr];

    for(len_t ir = 0; ir<nr; ir++){
        n1 = np1[ir] + (fluxGridType == FLUXGRIDTYPE_P1);
        n2 = np2[ir] + (fluxGridType == FLUXGRIDTYPE_P2);

        theta_bounceGrid[ir]   = new real_t*[n1*n2];
        weights_bounceGrid[ir] = new real_t*[n1*n2];
        B_bounceGrid[ir]       = new real_t*[n1*n2];
        ROverR0_bounceGrid[ir] = new real_t*[n1*n2];
        NablaR2_bounceGrid[ir] = new real_t*[n1*n2];

        for(len_t j = 0; j<n2; j++){
            for(len_t i = 0; i<n1; i++){
                pind = j*n1+i;
                theta_bounceGrid[ir][pind]   = new real_t[ntheta_interp_trapped];
                weights_bounceGrid[ir][pind] = new real_t[ntheta_interp_trapped];
                B_bounceGrid[ir][pind]       = new real_t[ntheta_interp_trapped];
                ROverR0_bounceGrid[ir][pind] = new real_t[ntheta_interp_trapped];
                NablaR2_bounceGrid[ir][pind] = new real_t[ntheta_interp_trapped];
            }
        }
    }


}

void BounceAverager::AllocateBounceIntegralQuantities(){
    DeallocateBounceIntegralQuantities();

    isTrapped    = new bool*[nr];
    isTrapped_fr = new bool*[nr+1];
    isTrapped_f1 = new bool*[nr];
    isTrapped_f2 = new bool*[nr];
    theta_b1     = new real_t*[nr];
    theta_b1_fr  = new real_t*[nr+1];
    theta_b1_f1  = new real_t*[nr];
    theta_b1_f2  = new real_t*[nr];
    theta_b2     = new real_t*[nr];
    theta_b2_fr  = new real_t*[nr+1];
    theta_b2_f1  = new real_t*[nr];
    theta_b2_f2  = new real_t*[nr];

    len_t n1, n2;
    for(len_t ir = 0; ir<nr; ir++){
        n1 = np1[ir];
        n2 = np2[ir];
        isTrapped[ir] = new bool[n1*n2];
        theta_b1[ir]  = new real_t[n1*n2];
        theta_b2[ir]  = new real_t[n1*n2];

        n1 = np1[ir]+1;
        n2 = np2[ir];    
        isTrapped_f1[ir] = new bool[n1*n2];
        theta_b1_f1[ir]  = new real_t[n1*n2];
        theta_b2_f1[ir]  = new real_t[n1*n2];

        n1 = np1[ir];
        n2 = np2[ir]+1;    
        isTrapped_f2[ir] = new bool[n1*n2];
        theta_b1_f2[ir]  = new real_t[n1*n2];
        theta_b2_f2[ir]  = new real_t[n1*n2];
    }
    for(len_t ir = 0; ir<=nr; ir++){
        n1 = np1[ir];
        n2 = np2[ir];
        isTrapped_fr[ir] = new bool[n1*n2];
        theta_b1_fr[ir]  = new real_t[n1*n2];
        theta_b2_fr[ir]  = new real_t[n1*n2];
    }
}


/**
 * Maybe we should hand these over to FVM::Grid? 
 */
void BounceAverager::DeallocateBounceIntegralQuantities(){
    if(isTrapped == nullptr)
        return;

    for(len_t ir = 0; ir<nr; ir++){
        delete [] isTrapped[ir];
        delete [] isTrapped_f1[ir];
        delete [] isTrapped_f2[ir];
        delete [] theta_b1[ir];
        delete [] theta_b1_f1[ir];
        delete [] theta_b1_f2[ir];
        delete [] theta_b2[ir];
        delete [] theta_b2_f1[ir];
        delete [] theta_b2_f2[ir];
    }
    for(len_t ir = 0; ir<=nr; ir++){
        delete [] isTrapped_fr[ir];
        delete [] theta_b1_fr[ir];
        delete [] theta_b2_fr[ir];
    }
    delete [] isTrapped;
    delete [] isTrapped_fr;
    delete [] isTrapped_f1;
    delete [] isTrapped_f2;
    delete [] theta_b1;
    delete [] theta_b1_fr;
    delete [] theta_b1_f1;
    delete [] theta_b1_f2;
    delete [] theta_b2;
    delete [] theta_b2_fr;
    delete [] theta_b2_f1;
    delete [] theta_b2_f2;
}


