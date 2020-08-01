/**
 * Implementation of BounceAverger class which handles everything 
 * needed to carry out bounce averages in DREAM. Each Grid object 
 * owns one instance of BounceAverager, which it uses to construct 
 * bounce averages on the grid.
 * 
 * Separate quadratures are used for trapped and passing orbits,
 * where the passing quadrature is inherited from the 
 * FluxSurfaceAverager provided in the constructor, and trapped
 * orbits use the.specified method (q_method_trapped). 
 * 
 * It is initialized (or refreshed) with Rebuild(), which should
 * be called after the FluxSurfaceAverager has been Rebuilt.
 */

#include "FVM/Grid/BounceAverager.hpp"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include <iostream>

using namespace std;
using namespace DREAM::FVM;

/**
 * Constructor.
 */
BounceAverager::BounceAverager(
    Grid *g, FluxSurfaceAverager* fsa, len_t ntheta_interp_trapped,
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

    // Use the Brent algorithm for root finding when determining the theta bounce points
    const gsl_root_fsolver_type *GSL_rootsolver_type = gsl_root_fsolver_brent;
    gsl_fsolver = gsl_root_fsolver_alloc (GSL_rootsolver_type);
    gsl_adaptive = gsl_integration_workspace_alloc(1000);
    gsl_acc = gsl_interp_accel_alloc();
    
}

/**
 * Destructor.
 */
BounceAverager::~BounceAverager(){
    gsl_interp_accel_free(gsl_acc);
    gsl_root_fsolver_free(gsl_fsolver);
    gsl_integration_workspace_free(gsl_adaptive);
    
    if(!integrateTrappedAdaptive)
        gsl_integration_fixed_free(gsl_w);
    delete [] np1;
    delete [] np2;
}


/**
 * Sets quantities related to the quadrature rule used for evaluating
 * bounce integrals on trapped orbits. We use a quadrature of the form
 *      integral( w(x, xmin, xmax)*F(x), xmin, xmax )
 *          = sum_i w_i F(x_i)
 * where theta is x_i, weights is w_i and quadratureWeightFunction 
 * is w(x, xmin, xmax).  ntheta_interp is the number of points i 
 * in the sum. The same quadrature is used on all trapped phase-space points.
 * Refer to GSL integration documentation for possible quadrature rules 
 * and corresponding weight functions, in case the options should be extended:
 *      https://www.gnu.org/software/gsl/doc/html/integration.html
 */
void BounceAverager::InitializeQuadrature(FluxSurfaceAverager::quadrature_method q_method){
    const gsl_integration_fixed_type *quadratureRule;
    std::function<real_t(real_t,real_t,real_t)> QuadFunc;
    switch(q_method){
        // Legendre quadrature is often best for smooth finite functions on a finite interval
        case FluxSurfaceAverager::QUAD_FIXED_LEGENDRE:
            quadratureRule = gsl_integration_fixed_legendre;
            QuadFunc = [](real_t /*x*/, real_t /*x_min*/, real_t /*x_max*/)
                                    {return 1;};
            break;
        /**
         * Chebyshev quadrature is best for integration along trapped orbits 
         * (unless the integrand is weighed by xi) where the metric takes the 
         * form of this weight function (square-root singularity on boundary).
         */
        case FluxSurfaceAverager::QUAD_FIXED_CHEBYSHEV:
            quadratureRule = gsl_integration_fixed_chebyshev;
            QuadFunc = [](real_t x, real_t x_min, real_t x_max)
                                    {return 1/sqrt((x_max-x)*(x-x_min) );};
            break;
        // Adaptive quadrature always works well but is significantly slower.
        case FluxSurfaceAverager::QUAD_ADAPTIVE:
            integrateTrappedAdaptive = true;
            return;
        default:
            throw FVMException("Quadrature rule '%d' not supported by BounceAverager.", q_method);            
    }
    /**
     * Generate reference quadrature on [0,1] which is later interpolated to the 
     * poloidal-angle interval [theta_b1, theta_b2] between the bounce points.
     */
    real_t unusedArgument = 0;
    real_t xmin = 0, xmax = 1;
    gsl_w = gsl_integration_fixed_alloc(quadratureRule,ntheta_interp_trapped,xmin,xmax,unusedArgument,unusedArgument);
    theta_trapped_ref = gsl_w->x;
    weights_trapped_ref = gsl_w->weights;
    // Divide by quadrature weight function to cast the integral on the desired form
    for(len_t it=0; it<ntheta_interp_trapped; it++)
        weights_trapped_ref[it] /= QuadFunc(theta_trapped_ref[it],0,1);

}


/**
 * Rebuilds quantities needed to perform bounce averages.
 * Should be called after FluxSurfaceAverager->Rebuild().
 */
void BounceAverager::Rebuild(){
    UpdateGridResolution();
    
    // If data has been generated already, deallocate.
    if(isTrapped!=nullptr){
        Metric->DeallocateData();
        B->DeallocateData();
        ROverR0->DeallocateData();
        NablaR2->DeallocateData();
    }

    // Initialize isTrapped and theta bounce points. True if any isTrapped.
    bool hasTrapped = InitializeBounceIntegralQuantities();

    // true if data should be stored on trapped grid (ie not adaptive quad) 
    bool storeTrapped =  hasTrapped && !integrateTrappedAdaptive;

    // true if data should be stored on passing grid 
    bool storePassing = !integratePassingAdaptive;

    // Allocate memory
    if ( storePassing || storeTrapped)
        Metric->AllocateData();
    if(storeTrapped){
        B->AllocateData();
        ROverR0->AllocateData();
        NablaR2->AllocateData();
    }

    // If using fixed quadrature on passing grid, store metric on passing theta grid
    // (the other quantities are already set via FluxSurfaceAverager)
    if(storePassing)
        Metric->SetDataForPassing(ntheta_interp_passing, theta_passing);
        
    // If using fixed quadrature on passing grid, store everything on trapped theta grid
    if (storeTrapped){
        B      ->SetDataForTrapped(ntheta_interp_trapped,theta_trapped_ref);
        ROverR0->SetDataForTrapped(ntheta_interp_trapped,theta_trapped_ref);
        NablaR2->SetDataForTrapped(ntheta_interp_trapped,theta_trapped_ref);
        Metric ->SetDataForTrapped(ntheta_interp_trapped,theta_trapped_ref);
    }

    // Calculate bounce-averaged metric and hand over to grid 
    real_t **Vp, **Vp_fr, **Vp_f1, **Vp_f2;
    SetVp(Vp,FLUXGRIDTYPE_DISTRIBUTION);
    SetVp(Vp_fr,FLUXGRIDTYPE_RADIAL);
    SetVp(Vp_f1,FLUXGRIDTYPE_P1);
    SetVp(Vp_f2,FLUXGRIDTYPE_P2);

    /**
     * TODO: This assumes that a PXi grid is used. It is probably harmless to create
     *  it if using ppar-pperp grid, in which case this quantity would never be used (?).       
     */
    real_t **VpOverP2AtZero = new real_t*[nr];
    for(len_t ir=0; ir<nr; ir++){
        VpOverP2AtZero[ir] = new real_t[np2[ir]];
        for(len_t j=0; j<np2[ir];j++){
            VpOverP2AtZero[ir][j] = grid->GetRadialGrid()->EvaluatePXiBounceIntegralAtP(ir,  0,  grid->GetMomentumGrid(ir)->GetP2(j),  FLUXGRIDTYPE_P1, [](real_t,real_t,real_t,real_t){return 1;});
        }
    }
    grid->SetVp(Vp,Vp_fr,Vp_f1,Vp_f2,VpOverP2AtZero);
}

/**
 *  Allocate and set VPrime to the bounce integral of the metric
 */
void BounceAverager::SetVp(real_t**&Vp, fluxGridType fluxGridType){
    function<real_t(real_t,real_t,real_t,real_t)> unityFunc 
               = [](real_t,real_t,real_t,real_t){return 1;};

    len_t nr = this->nr + (fluxGridType == FLUXGRIDTYPE_RADIAL);
    // XXX: assume same grid at all radii
    len_t n1 = np1[0] + (fluxGridType == FLUXGRIDTYPE_P1 ? 1 : 0);
    len_t n2 = np2[0] + (fluxGridType == FLUXGRIDTYPE_P2 ? 1 : 0);
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

/**
 * Evaluates the bounce average {F} of a function 
 *      F = F(xi/xi0, B/Bmin, R/R0, |nabla r|^2) on grid point (ir,i,j).
 */ 
real_t BounceAverager::CalculateBounceAverage(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F){
    real_t Vp = GetVp(ir,i,j,fluxGridType);    

    /**
     * Treat special case: 
     * Either r=0 (particle doesn't move in the poloidal plane) or 
     * xi0=0 in inhomogeneous magnetic field (infinitely deeply trapped on low-field side).
     * TODO: consider more carefully what happens with xi/xi0 in the xi0=0 case.
     */
    if(Vp == 0)  
        return F(1,1,1,1);

    // Otherwise, use regular definition:
    return EvaluateBounceIntegral(ir,i,j,fluxGridType, F) / Vp;
}


/**
 * The function returns the integrand of the bounce integral at theta.
 * Is used with adaptive quadrature.
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
    std::function<real_t(real_t,real_t,real_t,real_t)> F = params->Function;
    real_t Bmin = params->Bmin;
    real_t BOverBmin, xiOverXi0;
    if(B==Bmin){ 
        BOverBmin = 1;
        xiOverXi0 = 1;
    } else {
        BOverBmin = B/Bmin;
        real_t xi0Sq = xi0*xi0;
        xiOverXi0 = sqrt( (1 - BOverBmin*(1-xi0Sq))/xi0Sq );
    }
        
    return 2*M_PI*Metric*F(xiOverXi0,BOverBmin, ROverR0, NablaR2);
}

/**
 * Core function of this class: evaluates the bounce integral 
 *    BounceIntegral(X) = \int sqrt(g)*X dphi dtheta dzeta
 * taken over toroidal, poloidal and gyro angle, weighted by the 
 * phase-space metric sqrt(g). See doc/notes/theory section on 
 * Bounce averages for further details.
 */
real_t BounceAverager::EvaluateBounceIntegral(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F){
    real_t xi0 = GetXi0(ir,i,j,fluxGridType);
    real_t Bmin = fluxSurfaceAverager->GetBmin(ir,fluxGridType);
    
    bool isTrapped = BounceSurfaceQuantity::IsTrapped(ir,i,j,fluxGridType,grid);
    std::function<real_t(real_t,real_t,real_t,real_t)> F_eff;
    
    // If trapped, sum quantity over both directions along the field line
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
    real_t weightScaleFactor;
    if(isTrapped){
        ntheta = ntheta_interp_trapped;
        weights = weights_trapped_ref;
        weightScaleFactor = theta_b2-theta_b1;
    } else {
        ntheta = ntheta_interp_passing;
        weights = this->weights_passing;
        weightScaleFactor = 1;
    }

    real_t xiOverXi0,w,BOverBmin;        
    for (len_t it = 0; it<ntheta; it++) {
        // treat the singular cylindrical case 
        if(B[it]==Bmin){
            xiOverXi0 = 1;
            BOverBmin = 1;
        } else {
            real_t xi0Sq = xi0*xi0;
            BOverBmin = B[it]/Bmin;
            xiOverXi0 = sqrt((1- BOverBmin * (1-xi0Sq))/xi0Sq);
        }
        w = weightScaleFactor*weights[it];

        BounceIntegral += 2*M_PI*w*Metric[it]*F_eff(xiOverXi0,BOverBmin,ROverR0[it],NablaR2[it]);
    }        
    return BounceIntegral;
    
}



/**
 * Set isTrapped and poloidal bounce points and hands over ownership to Grid.
 */
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

/**
 * Helper function for InitializeBounceIntegralQuantities().
 */
bool BounceAverager::SetIsTrapped(bool **&isTrapped, real_t **&theta_b1, real_t **&theta_b2, fluxGridType fluxGridType){
    bool hasTrapped = false;
    real_t Bmin, Bmax;
    /**
     * XXX: Here we assume same grid at all radii.
     */
    len_t nr = this->nr + (fluxGridType == FLUXGRIDTYPE_RADIAL);
    len_t n1 = np1[0] + (fluxGridType == FLUXGRIDTYPE_P1);
    len_t n2 = np2[0]+ (fluxGridType == FLUXGRIDTYPE_P2);
    for(len_t ir = 0; ir<nr; ir++){
        real_t theta_Bmin = -1, theta_Bmax = -1;
        Bmin = fluxSurfaceAverager->GetBmin(ir,fluxGridType, &theta_Bmin);
        Bmax = fluxSurfaceAverager->GetBmax(ir,fluxGridType, &theta_Bmax);

        // in cylindrical grid or at r=0, isTrapped is false and we skip to next radius 
        if(Bmin==Bmax){
            for(len_t i = 0; i<n1*n2; i++)
                isTrapped[ir][i] = false;
            continue;
        }

        for(len_t j = 0; j<n2; j++)
            for(len_t i = 0; i<n1; i++){
                len_t pind = j*n1+i;
                real_t xi0 = GetXi0(ir,i,j,fluxGridType);
                if((1-xi0*xi0) > Bmin/Bmax){
                    isTrapped[ir][pind] = true;
                    hasTrapped = true;
                    FluxSurfaceAverager::FindBouncePoints(ir,Bmin, theta_Bmin, theta_Bmax, B->GetFluxSurfaceQuantity(), xi0, fluxGridType,
                                    &theta_b1[ir][pind],&theta_b2[ir][pind], gsl_fsolver);
                } else 
                    isTrapped[ir][n1*j+i] = false;
            }
    }
    return hasTrapped;
}


/**
 * Allocator.
 */
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

/**
 * Helper function to get Bmin from RadialGrid.
 */
/*real_t BounceAverager::GetBmin(len_t ir,fluxGridType fluxGridType){
    if (fluxGridType == FLUXGRIDTYPE_RADIAL){
        return grid->GetRadialGrid()->GetBmin_f(ir);
    } else {
        return grid->GetRadialGrid()->GetBmin(ir);
    }
}
*/
/**
 * Helper function to get Bmax from RadialGrid.
 */
/*
real_t BounceAverager::GetBmax(len_t ir,fluxGridType fluxGridType){
    if (fluxGridType == FLUXGRIDTYPE_RADIAL){
        return grid->GetRadialGrid()->GetBmax_f(ir);
    } else {
        return grid->GetRadialGrid()->GetBmax(ir);
    }
}
*/
/**
 * Helper function to get xi0 from MomentumGrid.
 */
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

/**
 * Helper function to get Vp from Grid.
 */
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

/**
 * Sets nr, np1 and np2.
 */
void BounceAverager::UpdateGridResolution(){
    nr = grid->GetNr();
    if(np1!=nullptr){
        delete [] np1;
        delete [] np2;
    }
    np1 = new len_t[nr];
    np2 = new len_t[nr];
    for(len_t ir = 0; ir<nr; ir++){
        np1[ir] = grid->GetNp1(ir);
        np2[ir] = grid->GetNp2(ir);
    }
    B->SetGridResolution(nr,np1,np2);
    ROverR0->SetGridResolution(nr,np1,np2);
    NablaR2->SetGridResolution(nr,np1,np2);
    Metric->SetGridResolution(nr,np1,np2);
}



/**
 * Helper function for avalanche deltaHat calculation: 
 * returns the function xi0Star = sqrt( 1 - 1/BOverBmin * 2/(g+1) )
 * See documentation in doc/notes/theory.
 */
real_t xi0Star(real_t BOverBmin, real_t p){
    real_t g = sqrt(1+p*p);
    // This form is numerically stable for arbitrary p and BOverBmin
    return sqrt( (p*p/(g+1) + 2*(BOverBmin-1)/BOverBmin ) / (g+1) );
}

/**
 * For gsl root finding: returns the function sgn*(xi0Star - xi0), which 
 * defines the integration limits in avalanche deltaHat calculation.
 */
struct xiStarParams {real_t p; real_t xi0; len_t ir; real_t Bmin; const BounceSurfaceQuantity *B; real_t sgn;};
real_t xi0StarRootFunc(real_t theta, void *par){
    struct xiStarParams *params = (struct xiStarParams *) par;
    len_t ir = params->ir;
    const BounceSurfaceQuantity *B = params->B;
    real_t Bmin = params->Bmin;
    real_t sgn = params->sgn;
    real_t BOverBmin = B->evaluateAtTheta(ir, theta, FLUXGRIDTYPE_DISTRIBUTION)/Bmin;
    real_t p = params->p;
    real_t xi0 = params->xi0;
    return sgn*(xi0Star(BOverBmin,p) - xi0);
}


struct hParams {real_t p; len_t ir; real_t Bmin; real_t Vp; real_t dxi; const BounceSurfaceQuantity *B;};
real_t hIntegrand(real_t theta, void *par){
    struct hParams *params = (struct hParams *) par;
    len_t ir = params->ir;
    const BounceSurfaceQuantity *B_qty = params->B;
    real_t Bmin = params->Bmin;
    real_t B = B_qty->evaluateAtTheta(ir, theta, FLUXGRIDTYPE_DISTRIBUTION);
    real_t BOverBmin = B/Bmin;
    real_t p = params->p;
    real_t Vp = params->Vp;
    real_t dxi = params->dxi;

    real_t g = sqrt(1+p*p);
    real_t xi = sqrt((g-1)/(g+1));
    real_t xi0 = xi0Star(BOverBmin,p);
    real_t sqrtgOverP2 = MomentumGrid::evaluatePXiMetricOverP2(p,xi0,B,Bmin);
    // 2*pi for the trivial phi integral
    return 2*M_PI * xi/xi0 * sqrtgOverP2 / (dxi * Vp);
}


real_t BounceAverager::EvaluateAvalancheDeltaHat(len_t ir, real_t p, real_t xi_l, real_t xi_u, real_t Vp, real_t VpVol){
    real_t theta_Bmin=0, theta_Bmax=0;
    real_t Bmin = fluxSurfaceAverager->GetBmin(ir, FLUXGRIDTYPE_DISTRIBUTION,&theta_Bmin);
    real_t Bmax = fluxSurfaceAverager->GetBmax(ir, FLUXGRIDTYPE_DISTRIBUTION,&theta_Bmax);
    if( xi0Star(Bmax/Bmin,p) <= xi_l )
        return 0;
    else if( xi0Star(1,p) >= xi_u )
        return 0;

    // Return something reasonable for Vp=0 (corresponding to infinitely deeply trapped particles xi0=0... 
    // this point is singular and there is no clearcut way of defining it -- 
    // in the end they shouldn't contribute anyway)
    if(Vp==0)
        return 0; //placeholder
//        return 1 / (p*p*(xi_u-xi_l));

    // PLACEHOLDER: assumes that the source is a delta function also in xi0,
    // which approximates the pitch distribution of the secondaries. General
    // calculation below is TODO
    return 2*M_PI*VpVol/(Vp/(p*p)) *  grid->GetRadialGrid()->GetFSA_B(ir) / (p*p*(xi_u-xi_l));


    // else, there are two nontrivial intervals [theta_l, theta_u] on which contributions are obtained


    xiStarParams xi_params_u = {p,xi_u,ir,Bmin,B, -1.0}; 
    xiStarParams xi_params_l = {p,xi_l,ir,Bmin,B, 1.0}; 

    bool upperForAllTheta = (xi0StarRootFunc(theta_Bmin, &xi_params_u) > 0) && (xi0StarRootFunc(theta_Bmax, &xi_params_u) > 0);
    bool lowerForAllTheta = (xi0StarRootFunc(theta_Bmin, &xi_params_l) > 0) && (xi0StarRootFunc(theta_Bmax, &xi_params_l) > 0);

    // if all poloidal angles contribute fully to the integral, return the known exact value.
    if(upperForAllTheta && lowerForAllTheta)
        return VpVol/(Vp/(p*p))*grid->GetRadialGrid()->GetFSA_B(ir) / (p*p*(xi_u-xi_l));


    gsl_function gsl_func;
        

    /**
    * theta2 is assumed to exist on the interval [theta_Bmin, theta_Bmax], 
    * and theta1 on the interval [theta_Bmax-2*pi, theta_Bmin]
    */
    real_t theta_u1, theta_u2, theta_l1, theta_l2;
    if(upperForAllTheta){
        theta_u1 = theta_Bmax - 2*M_PI;
        theta_u2 = theta_Bmin;
    } else {
        gsl_func.function = &(xi0StarRootFunc);
        gsl_func.params = &xi_params_u;
        FluxSurfaceAverager::FindThetas(theta_Bmin,theta_Bmax,&theta_u1, &theta_u2, gsl_func, gsl_fsolver);
    }
    
    if(lowerForAllTheta){
        theta_l1 = theta_Bmax;
        theta_l2 = theta_Bmin;
    } else {
        gsl_func.function = &(xi0StarRootFunc);
        gsl_func.params = &xi_params_l;
        FluxSurfaceAverager::FindThetas(theta_Bmin,theta_Bmax,&theta_l1, &theta_l2, gsl_func, gsl_fsolver);
    }


    // if the midpoint between theta_l1 and theta_u1 does not satisfy the condition,
    // theta_u1 and theta_u2 have been mixed up in the root finding algorithm and 
    // should be switched.
    bool midpointIsValid = (xi0StarRootFunc( (theta_l1+theta_u1)/2, &xi_params_u ) > 0 );
    if(!midpointIsValid){
        real_t tmp = theta_u1;
        theta_u1 = theta_u2;
        theta_u2 = tmp;
    }

//struct hParams {real_t p; len_t ir; real_t Bmin; real_t Vp; real_t dxi; FluxSurfaceAverager *fsa; const BounceSurfaceQuantity *B;};

    hParams h_params = {p,ir,Bmin,Vp,xi_u-xi_l, B};
    gsl_func.function = &(hIntegrand);
    gsl_func.params = &h_params;
    
    real_t deltaHat1, deltaHat2;
    real_t epsabs = 0, epsrel = 1e-4, lim = gsl_adaptive->limit, error;
    gsl_integration_qags(&gsl_func, min(theta_l1,theta_u1), max(theta_l1,theta_u1),epsabs,epsrel,lim,gsl_adaptive,&deltaHat1, &error);
    gsl_integration_qags(&gsl_func, min(theta_l2,theta_u2), max(theta_l2,theta_u2),epsabs,epsrel,lim,gsl_adaptive,&deltaHat2, &error);



    // integrate h on [theta_l1, theta_u1] and [theta_l2, theta_u2]

    // placeholder: cylindrical limit also in the inhomogeneous case
    return deltaHat1+deltaHat2;
}
