/**
 * Implementation of the FluxSurfaceAverager class which 
 * handles everything needed to carry out flux surface 
 * averages in DREAM. Is initialized by a RadialGridGenerator
 * via the method SetReferenceMagneticFieldData.
 */

#include "FVM/Grid/BounceAverager.hpp"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

using namespace std;
using namespace DREAM::FVM;

/**
 * Constructor.
 */
BounceAverager::BounceAverager(
    Grid *g, FluxSurfaceAverager* fsa, len_t ntheta_interp_trapped,
//    enum OptionConstants::momentumgrid_type mgtype,
    FluxSurfaceAverager::quadrature_method q_method_trapped
) : grid(g), fluxSurfaceAverager(fsa), ntheta_interp_trapped(ntheta_interp_trapped) {

    geometryIsSymmetric = fluxSurfaceAverager->isGeometrySymmetric();
    ntheta_interp_passing       = fluxSurfaceAverager->GetNTheta();
    theta_passing               = fluxSurfaceAverager->GetTheta();
    weights_passing             = fluxSurfaceAverager->GetWeights();
    integratePassingAdaptive    = fluxSurfaceAverager->isIntegrationAdaptive();

    B       = new BounceSurfaceQuantity(grid, fluxSurfaceAverager->GetB());
    ROverR0 = new BounceSurfaceQuantity(grid, fluxSurfaceAverager->GetROverR0());
    NablaR2 = new BounceSurfaceQuantity(grid, fluxSurfaceAverager->GetNablaR2());
    Metric  = new BounceSurfaceMetric(grid, fluxSurfaceAverager->GetJacobian(), fluxSurfaceAverager->GetB());
    

    InitializeQuadrature(q_method_trapped);

    // Use the Brent algorithm for root finding in determining the theta bounce points
    const gsl_root_fsolver_type *GSL_rootsolver_type = gsl_root_fsolver_brent;
    gsl_fsolver = gsl_root_fsolver_alloc (GSL_rootsolver_type);

    gsl_acc = gsl_interp_accel_alloc();
    
}

BounceAverager::~BounceAverager(){
    gsl_interp_accel_free(gsl_acc);
    gsl_root_fsolver_free(gsl_fsolver);
    
    if(integrateTrappedAdaptive)
        gsl_integration_workspace_free(gsl_adaptive);
    else
        gsl_integration_fixed_free(gsl_w);
    delete [] np1;
    delete [] np2;
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
            break;
        case FluxSurfaceAverager::QUAD_FIXED_CHEBYSHEV:
            // Chebyshev quadrature may be better for integration along trapped orbits
            // where the metric takes the form of this weight function.
            quadratureRule = gsl_integration_fixed_chebyshev;
            QuadFunc = [](real_t x, real_t x_min, real_t x_max)
                                    {return 1/sqrt((x_max-x)*(x-x_min) );};
            break;
        case FluxSurfaceAverager::QUAD_ADAPTIVE:
            gsl_adaptive = gsl_integration_workspace_alloc(1000);
            integrateTrappedAdaptive = true;
            return;
        default:
            throw FVMException("Quadrature rule '%d' not supported by FluxSurfaceAverager.", q_method);            
    }
    real_t unusedArgument = 0;
    real_t xmin = 0, xmax = 1;
    gsl_w = gsl_integration_fixed_alloc(quadratureRule,ntheta_interp_trapped,xmin,xmax,unusedArgument,unusedArgument);
    theta_trapped_ref = gsl_w->x;
    weights_trapped_ref = gsl_w->weights;
    for(len_t it=0; it<ntheta_interp_trapped; it++)
        weights_trapped_ref[it] /= QuadFunc(theta_trapped_ref[it],0,1);

}


/**
 * Rebuilds quantities needed to perform bounce averages.
 * Should be called after FluxSurfaceAverager->Rebuild().
 * Calculates quantities which are independent of poloidal
 * angle theta, and hands them over to Grid.
 */
void BounceAverager::Rebuild(){
    UpdateGridResolution();
    
    if(isTrapped==nullptr){
        Metric->DeallocateData();
        B->DeallocateData();
        ROverR0->DeallocateData();
        NablaR2->DeallocateData();
    }

    bool hasTrapped = InitializeBounceIntegralQuantities();

    // true if data should be stored on fixed-quadrature trapped grid 
    bool storeTrapped =  hasTrapped && !integrateTrappedAdaptive;

    // true if data should be stored on fixed-quadrature passing grid 
    bool storePassing = !integratePassingAdaptive;

    // allocate metric data 
    if ( storePassing || storeTrapped)
        Metric->AllocateData();
    if(storeTrapped){
        B->AllocateData();
        ROverR0->AllocateData();
        NablaR2->AllocateData();
    }

    // If using fixed quadrature on passing grid, store metric on passing theta grid
    if(storePassing)
        Metric->SetDataForPassing(ntheta_interp_passing, theta_passing);
        
    // If using fixed quadrature on passing grid, store everything on trapped theta grid
    if (storeTrapped){
        B      ->SetDataForTrapped(ntheta_interp_trapped,theta_trapped_ref);
        ROverR0->SetDataForTrapped(ntheta_interp_trapped,theta_trapped_ref);
        NablaR2->SetDataForTrapped(ntheta_interp_trapped,theta_trapped_ref);
        Metric ->SetDataForTrapped(ntheta_interp_trapped,theta_trapped_ref);
    }


    real_t **Vp, **Vp_fr, **Vp_f1, **Vp_f2;
    SetVp(Vp,FLUXGRIDTYPE_DISTRIBUTION);
    SetVp(Vp_fr,FLUXGRIDTYPE_RADIAL);
    SetVp(Vp_f1,FLUXGRIDTYPE_P1);
    SetVp(Vp_f2,FLUXGRIDTYPE_P2);

    real_t **VpOverP2AtZero = new real_t*[nr];
    for(len_t ir=0; ir<nr; ir++){
        VpOverP2AtZero[ir] = new real_t[np2[ir]];
        for(len_t j=0; j<np2[ir];j++){
            VpOverP2AtZero[ir][j] = grid->GetRadialGrid()->EvaluatePXiBounceIntegralAtP(ir,  0,  grid->GetMomentumGrid(ir)->GetP2(j),  FLUXGRIDTYPE_P1, [](real_t,real_t,real_t,real_t){return 1;});
        }
    }
    grid->SetVp(Vp,Vp_fr,Vp_f1,Vp_f2,VpOverP2AtZero);
}

void BounceAverager::SetVp(real_t**&Vp, fluxGridType fluxGridType){
    function<real_t(real_t,real_t,real_t,real_t)> unityFunc 
               = [](real_t,real_t,real_t,real_t){return 1;};

    len_t nr = this->nr + (fluxGridType == FLUXGRIDTYPE_RADIAL);
    // XXX: assume same grid at all radii
    len_t n1 = np1[0] + (fluxGridType == FLUXGRIDTYPE_P1);
    len_t n2 = np2[0] + (fluxGridType == FLUXGRIDTYPE_P2);
    Vp = new real_t*[nr];
    for(len_t ir = 0; ir<nr; ir++){
        Vp[ir] = new real_t[n1*n2];
        for(len_t j = 0; j<n2; j++)
            for(len_t i = 0; i<n1; i++){
                len_t pind = j*n1+i;
                Vp[ir][pind] = EvaluateBounceIntegral(ir,i,j,fluxGridType,unityFunc);
            }
    }
}

// Evaluates the bounce average {F} of a function F = F(xi/xi0, B/Bmin, R/R0, |nabla r|^2) on grid point (ir,i,j). 
real_t BounceAverager::CalculateBounceAverage(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F){
    /**
     * XXX: Here we assume same grid at all radii.
     */
    real_t xi0 = GetXi0(ir,i,j,fluxGridType);
    real_t Vp = GetVp(ir,i,j,fluxGridType);    
    
    if(xi0*xi0 < 1e-30){
        real_t BounceAverageSingular, sqrtterm;
        for(len_t it=0; it<ntheta_interp_trapped; it++){
            sqrtterm = sqrt(1-theta_trapped_ref[it]*theta_trapped_ref[it]);
            BounceAverageSingular +=  weights_trapped_ref[it]* (F(sqrtterm,1,1,1)+F(-sqrtterm,1,1,1))/sqrtterm;
        }
        BounceAverageSingular /= M_PI;
        return BounceAverageSingular;
    } else
        return EvaluateBounceIntegral(ir,i,j,fluxGridType, F) / Vp;
}


/**
 * The function returns the integrand of the bounce integral at theta
 */
struct BounceIntegralParams {
    real_t xi0; std::function<real_t(real_t,real_t,real_t,real_t)> Function; len_t ir; len_t i; len_t j;
    real_t Bmin; fluxGridType fgType; BounceAverager *bAvg;
};
real_t BounceAverager::BounceIntegralFunction(real_t theta, void *p){
    struct BounceIntegralParams *params = (struct BounceIntegralParams *) p;
    len_t ir = params->ir;
    len_t i = params->i;
    len_t j = params->j;
    real_t xi0 = params->xi0;
    fluxGridType fluxGridType = params->fgType;
    BounceAverager *bounceAverager = params->bAvg;
    real_t B = bounceAverager->GetB()->evaluateAtTheta(ir, theta, fluxGridType);
    real_t Metric  = bounceAverager->GetMetric()->evaluateAtTheta(ir, i, j, theta, fluxGridType);
    real_t ROverR0 = bounceAverager->GetROverR0()->evaluateAtTheta(ir, theta, fluxGridType);
    real_t NablaR2 = bounceAverager->GetNablaR2()->evaluateAtTheta(ir, theta, fluxGridType);
    real_t BOverBmin = B/params->Bmin;
    std::function<real_t(real_t,real_t,real_t,real_t)> F = params->Function;
    real_t xi0Sq = xi0*xi0;
    real_t xiOverXi0 = sqrt( (1 - BOverBmin*(1-xi0Sq))/xi0Sq );
    return 2*M_PI*Metric*F(xiOverXi0,BOverBmin, ROverR0, NablaR2);
}


real_t BounceAverager::EvaluateBounceIntegral(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F){
    real_t Bmin = GetBmin(ir,fluxGridType);
    real_t Bmax = GetBmax(ir,fluxGridType);
    bool BminEqBmax = (Bmin==Bmax);
    real_t xi0 = GetXi0(ir,i,j,fluxGridType);

    bool isTrapped = BounceSurfaceQuantity::IsTrapped(ir,i,j,fluxGridType,grid);
    std::function<real_t(real_t,real_t,real_t,real_t)> F_eff;
    // If trapped, adds contribution from -xi0, since negative xi0 are presumably not kept on the grid.
    if (isTrapped)
        F_eff = [&](real_t x, real_t  y, real_t z, real_t w){return  F(x,y,z,w) + F(-x,y,z,w) ;};
    else 
        F_eff = F;

    real_t theta_b1 = BounceSurfaceQuantity::Theta_B1(ir,i,j,fluxGridType,grid);
    real_t theta_b2 = BounceSurfaceQuantity::Theta_B2(ir,i,j,fluxGridType,grid);

    real_t BounceIntegral = 0;
    // If using adaptive-integration setting, perform bounce integral with GSL qags
    if( ( (!isTrapped) && (integratePassingAdaptive)) || (isTrapped && integrateTrappedAdaptive) ){
        gsl_function GSL_func; 
        BounceIntegralParams params = {xi0, F, ir, i, j, Bmin, fluxGridType, this}; 
        GSL_func.function = &(BounceIntegralFunction);
        GSL_func.params = &params;
        real_t epsabs = 0, epsrel = 1e-4, lim = gsl_adaptive->limit, error;
//        len_t key = GSL_INTEG_GAUSS31; // intermediate-order method in qag
        gsl_integration_qags(&GSL_func, theta_b1, theta_b2,epsabs,epsrel,lim,gsl_adaptive,&BounceIntegral, &error);
        return BounceIntegral;
    }
    // otherwise continue and use the chosen fixed quadrature

    const real_t *B       = this->B->GetData(ir,i,j,fluxGridType);
    const real_t *ROverR0 = this->ROverR0->GetData(ir,i,j,fluxGridType);
    const real_t *NablaR2 = this->NablaR2->GetData(ir,i,j,fluxGridType);
    const real_t *Metric  = this->Metric->GetData(ir,i,j,fluxGridType);

    len_t ntheta;
    const real_t *weights;
    if(isTrapped){
        ntheta = ntheta_interp_trapped;
        weights = weights_trapped_ref;
    } else {
        ntheta = ntheta_interp_passing;
        weights = this->weights_passing;
        theta_b1 = 0;
        theta_b2 = 1; // do not rescale weights in sum below
    }

    real_t xiOverXi0,w,BOverBmin;        
    for (len_t it = 0; it<ntheta; it++) {
        if(BminEqBmax){
            xiOverXi0 = 1;
            BOverBmin = 1;
        } else {
            real_t xi0Sq = xi0*xi0;
            xiOverXi0 = sqrt((1- B[it]/Bmin * (1-xi0Sq))/xi0Sq);
            BOverBmin = B[it]/Bmin;
        }
        w = (theta_b2-theta_b1)*weights[it];

        BounceIntegral += 2*M_PI*w*Metric[it]*F_eff(xiOverXi0,BOverBmin,ROverR0[it],NablaR2[it]);
    }        
    return BounceIntegral;
    
}



// Sets isTrapped and poloidal bounce points and hands over ownership to Grid.
bool BounceAverager::InitializeBounceIntegralQuantities(){
    AllocateBounceIntegralQuantities();
    bool hasTrapped;
    hasTrapped  = SetIsTrapped(isTrapped,    theta_b1,    theta_b2,    FLUXGRIDTYPE_DISTRIBUTION);
    hasTrapped += SetIsTrapped(isTrapped_fr, theta_b1_fr, theta_b2_fr, FLUXGRIDTYPE_RADIAL);
    hasTrapped += SetIsTrapped(isTrapped_f1, theta_b1_f1, theta_b2_f1, FLUXGRIDTYPE_P1);
    hasTrapped += SetIsTrapped(isTrapped_f2, theta_b1_f2, theta_b2_f2, FLUXGRIDTYPE_P2);

    grid->SetBounceParameters(
        isTrapped, isTrapped_fr, isTrapped_f1, isTrapped_f2,
        theta_b1,  theta_b1_fr,  theta_b1_f1,  theta_b1_f2, 
        theta_b2,  theta_b2_fr,  theta_b2_f1,  theta_b2_f2);

    return hasTrapped;
}


bool BounceAverager::SetIsTrapped(bool **&isTrapped, real_t **&theta_b1, real_t **&theta_b2, fluxGridType fluxGridType){
    bool hasTrapped = false;
    real_t Bmin, Bmax;
    len_t nr = this->nr + (fluxGridType == FLUXGRIDTYPE_RADIAL);
    len_t n1, n2;
    for(len_t ir = 0; ir<nr; ir++){
        /**
         * XXX: Here we assume same grid at all radii.
         */
        RadialGrid *rGrid = grid->GetRadialGrid();
        if (fluxGridType == FLUXGRIDTYPE_RADIAL){
            Bmin = rGrid->GetBmin_f(ir);
            Bmax = rGrid->GetBmax_f(ir);
        } else {
            Bmin = rGrid->GetBmin(ir);
            Bmax = rGrid->GetBmax(ir);
        }

        n1 = np1[0] + (fluxGridType == FLUXGRIDTYPE_P1);
        n2 = np2[0]+ (fluxGridType == FLUXGRIDTYPE_P2);
        // in cylindrical grid or at r=0, isTrapped is false and we skip to next radius 
        if(Bmin==Bmax){
            for(len_t i = 0; i<n1*n2; i++)
                isTrapped[ir][i] = false;
            continue;
        }

        for(len_t j = 0; j<n2; j++){
            for(len_t i = 0; i<n1; i++){
                len_t pind = j*n1+i;
                real_t xi0 = GetXi0(ir,i,j,fluxGridType);
                if((1-xi0*xi0) > Bmin/Bmax){
                    isTrapped[ir][pind] = true;
                    hasTrapped = true;
                    FluxSurfaceAverager::FindBouncePoints(ir,Bmin, B->GetFluxSurfaceQuantity(), xi0, fluxGridType,
                                    &theta_b1[ir][pind],&theta_b2[ir][pind], gsl_fsolver, geometryIsSymmetric);
                } else {
                    isTrapped[ir][n1*j+i] = false;

                }
            }
        }
    }
    return hasTrapped;
}



void BounceAverager::AllocateBounceIntegralQuantities(){
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
    // XXX: Assume same momentum grid at all radii
    for(len_t ir = 0; ir<=nr; ir++){
        n1 = np1[0];
        n2 = np2[0];
        isTrapped_fr[ir] = new bool[n1*n2];
        theta_b1_fr[ir]  = new real_t[n1*n2];
        theta_b2_fr[ir]  = new real_t[n1*n2];
    }
}


real_t BounceAverager::GetBmin(len_t ir,fluxGridType fluxGridType){
    if (fluxGridType == FLUXGRIDTYPE_RADIAL){
        return grid->GetRadialGrid()->GetBmin_f(ir);
    } else {
        return grid->GetRadialGrid()->GetBmin(ir);
    }
}
real_t BounceAverager::GetBmax(len_t ir,fluxGridType fluxGridType){
    if (fluxGridType == FLUXGRIDTYPE_RADIAL){
        return grid->GetRadialGrid()->GetBmax_f(ir);
    } else {
        return grid->GetRadialGrid()->GetBmax(ir);
    }
}

real_t BounceAverager::GetXi0(len_t ir, len_t i, len_t j, fluxGridType fluxGridType)
{
    if (fluxGridType == FLUXGRIDTYPE_P1) 
        return grid->GetMomentumGrid(ir)->GetXi0_f1(i,j);
    else if (fluxGridType == FLUXGRIDTYPE_P2) 
        return grid->GetMomentumGrid(ir)->GetXi0_f2(i,j);
    else if (fluxGridType == FLUXGRIDTYPE_RADIAL) // XXX: Assume same momentum grid at all radii
        return grid->GetMomentumGrid(0)->GetXi0(i,j);
    else
        return grid->GetMomentumGrid(ir)->GetXi0(i,j);        
}

real_t BounceAverager::GetVp(len_t ir, len_t i, len_t j, fluxGridType fluxGridType)
{
    if (fluxGridType == FLUXGRIDTYPE_P1) {
        return grid->GetVp_f1(ir,i,j);
    } else if (fluxGridType == FLUXGRIDTYPE_P2) {
        return grid->GetVp_f2(ir,i,j);
    } else if (fluxGridType == FLUXGRIDTYPE_RADIAL) { 
        return grid->GetVp_fr(ir,i,j);
    } else{
        return grid->GetVp(ir,i,j);
    }
}


void BounceAverager::UpdateGridResolution(){
    nr = grid->GetNr();
    if(np1!=nullptr){
        delete [] np1;
        delete [] np2;
    }
    np1 = new len_t[nr];
    np2 = new len_t[nr];
    for(len_t ir = 0; ir<nr; ir++){
        MomentumGrid *mg = grid->GetMomentumGrid(ir);
        np1[ir] = mg->GetNp1();
        np2[ir] = mg->GetNp2();
    }
    B->SetGridResolution(nr,np1,np2);
    ROverR0->SetGridResolution(nr,np1,np2);
    NablaR2->SetGridResolution(nr,np1,np2);
    Metric->SetGridResolution(nr,np1,np2);
}