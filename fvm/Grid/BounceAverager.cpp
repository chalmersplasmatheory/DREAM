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

    BOverBmin = new BounceSurfaceQuantity(grid, fluxSurfaceAverager->GetBOverBmin());
    ROverR0   = new BounceSurfaceQuantity(grid, fluxSurfaceAverager->GetROverR0());
    NablaR2   = new BounceSurfaceQuantity(grid, fluxSurfaceAverager->GetNablaR2());
    Metric    = new BounceSurfaceMetric(grid, fluxSurfaceAverager->GetJacobian(), fluxSurfaceAverager->GetBOverBmin(),fluxSurfaceAverager);

    InitializeQuadrature(q_method_trapped);

    // Use the Brent algorithm for root finding when determining the theta bounce points
    const gsl_root_fsolver_type *GSL_rootsolver_type = gsl_root_fsolver_brent;
    gsl_fsolver = gsl_root_fsolver_alloc (GSL_rootsolver_type);
    gsl_adaptive = gsl_integration_workspace_alloc(1000);
    gsl_acc = gsl_interp_accel_alloc();
    qaws_table = gsl_integration_qaws_table_alloc(-0.5, -0.5, 0, 0);
    
}

/**
 * Destructor.
 */
BounceAverager::~BounceAverager(){
    gsl_interp_accel_free(gsl_acc);
    gsl_root_fsolver_free(gsl_fsolver);
    gsl_integration_workspace_free(gsl_adaptive);
    gsl_integration_qaws_table_free(qaws_table);
    
    if(!integrateTrappedAdaptive)
        gsl_integration_fixed_free(gsl_w);
    delete [] np1;
    delete [] np2;
}

/** 
 * Calculations related to the bounce average of the delta function
 * appearing in the Rosenbluth-Putvinski avalanche source
 */
#include "BounceAverager.avalancheDelta.cpp"

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
    function<real_t(real_t,real_t,real_t)> QuadFunc;
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
        // Adaptive quadrature always works well but can be significantly slower.
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
        BOverBmin->DeallocateData();
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
        BOverBmin->AllocateData();
        ROverR0->AllocateData();
        NablaR2->AllocateData();
    }

    // If using fixed quadrature on passing grid, store metric on passing theta grid
    // (the other quantities are already set via FluxSurfaceAverager)
    if(storePassing)
        Metric->SetDataForPassing(ntheta_interp_passing, theta_passing);
        
    // If using fixed quadrature on passing grid, store everything on trapped theta grid
    if (storeTrapped){
        BOverBmin->SetDataForTrapped(ntheta_interp_trapped,theta_trapped_ref);
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
    function<real_t(real_t,real_t,real_t,real_t)> unityFunc =  [](real_t,real_t,real_t,real_t){return 1;};
    int_t unityList[5] = {0,0,0,0,1};
    real_t **VpOverP2AtZero = new real_t*[nr];
    for(len_t ir=0; ir<nr; ir++){
        VpOverP2AtZero[ir] = new real_t[np2[ir]];
        for(len_t j=0; j<np2[ir];j++)
            VpOverP2AtZero[ir][j] = grid->GetRadialGrid()->EvaluatePXiBounceIntegralAtP(
                ir, grid->GetMomentumGrid(ir)->GetP2(j),  
                FLUXGRIDTYPE_P1, unityFunc, unityList
            );
    }
    grid->SetVp(Vp,Vp_fr,Vp_f1,Vp_f2,VpOverP2AtZero);
}

/**
 *  Allocate and set VPrime to the bounce integral of the metric
 */
void BounceAverager::SetVp(real_t**&Vp, fluxGridType fluxGridType){
    function<real_t(real_t,real_t,real_t,real_t)> unityFunc 
               = [](real_t,real_t,real_t,real_t){return 1;};
    int_t unityList[5] = {0,0,0,0,1};
    len_t nr = this->nr + (fluxGridType == FLUXGRIDTYPE_RADIAL);
    // XXX: assume same grid at all radii
    len_t n1 = np1[0] + (fluxGridType == FLUXGRIDTYPE_P1 ? 1 : 0);
    len_t n2 = np2[0] + (fluxGridType == FLUXGRIDTYPE_P2 ? 1 : 0);
    Vp = new real_t*[nr];
    // XXX: assumes p-xi grid
    const real_t *p;
    if(fluxGridType == FLUXGRIDTYPE_P1)
        p = grid->GetMomentumGrid(0)->GetP1_f();
    else
        p = grid->GetMomentumGrid(0)->GetP1();

    bool isPXiGrid = true; 
    for(len_t ir = 0; ir<nr; ir++){
        Vp[ir] = new real_t[n1*n2];
        for(len_t j = 0; j<n2; j++)
            if(isPXiGrid){ // assume Vp scales as ~p^2
                len_t iOffset = (fluxGridType == FLUXGRIDTYPE_P1 ? 1 : 0);
                real_t Vp0 = EvaluateBounceIntegral(ir,iOffset,j,fluxGridType,unityFunc, unityList);
                real_t p0Sq = p[iOffset]*p[iOffset];
                for(len_t i = 0; i<n1; i++)
                    Vp[ir][j*n1+i] = Vp0 * p[i]*p[i] / p0Sq;
            } else 
                for(len_t i = 0; i<n1; i++)
                    Vp[ir][j*n1+i] = EvaluateBounceIntegral(ir,i,j,fluxGridType,unityFunc, unityList);
    }
}

/**
 * Evaluates the bounce average {F} of a function 
 *      F = F(xi/xi0, B/Bmin, R/R0, |nabla r|^2) on grid point (ir,i,j).
 */ 
real_t BounceAverager::CalculateBounceAverage(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, function<real_t(real_t,real_t,real_t,real_t)> F, int_t *F_list){
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
    return EvaluateBounceIntegral(ir,i,j,fluxGridType, F, F_list) / Vp;
}


/**
 * The function returns the integrand of the bounce integral at theta.
 * Is used with adaptive quadrature.
 */
struct BounceIntegralParams {
    real_t xi0; real_t p; function<real_t(real_t,real_t,real_t,real_t)> Function; int_t *F_list; len_t ir; 
    real_t theta_b1; real_t theta_b2; real_t Bmin; fluxGridType fgType; FluxSurfaceAverager *FSA; bool integrateQAWS;
};
real_t BounceAverager::BounceIntegralFunction(real_t theta, void *par){
    struct BounceIntegralParams *params = (struct BounceIntegralParams *) par;
    len_t ir = params->ir;

    real_t theta_b1 = params->theta_b1;
    real_t theta_b2 = params->theta_b2;
    real_t xi0 = params->xi0;
    real_t p = params->p;

    fluxGridType fluxGridType = params->fgType;
    FluxSurfaceAverager *FSA = params->FSA;

    real_t B,Jacobian,ROverR0,NablaR2;
    FSA->GeometricQuantitiesAtTheta(ir,theta,B,Jacobian,ROverR0,NablaR2,fluxGridType);
    real_t Bmin = params->Bmin;
    real_t BOverBmin=1;
    if(Bmin)
        BOverBmin = B/Bmin;
    real_t Metric = p*p*Jacobian*MomentumGrid::evaluatePXiMetricOverP2(xi0, BOverBmin);

    real_t xi2 = 1-BOverBmin*(1-xi0*xi0);
    real_t xiOverXi0 = 1;
    if(xi2>0)
        xiOverXi0 = sqrt(xi2/(xi0*xi0));

    int_t *Flist = params->F_list;
    real_t Function;
    if(Flist != nullptr) // if the function exponent list is provided, prioritize to build function this way
        Function = FSA->AssembleBAFunc(xiOverXi0, BOverBmin, ROverR0, NablaR2,Flist);
    else 
        Function = params->Function(xiOverXi0,BOverBmin,ROverR0,NablaR2);

    real_t S = 2*M_PI*Metric*Function;
    if(params->integrateQAWS) // divide by weight function in QAWS quadrature
        return S*sqrt((theta-theta_b1)*(theta_b2-theta));
    else 
        return S;
}


/**
 * Core function of this class: evaluates the bounce integral 
 *    BounceIntegral(X) = \int sqrt(g)*X dphi dtheta dzeta
 * taken over toroidal, poloidal and gyro angle, weighted by the 
 * phase-space metric sqrt(g). See doc/notes/theory section on 
 * Bounce averages for further details.
 */
const real_t NEAR_TRAPPED_BOUNDARY_THRESHOLD = 0.01; // distance in xi from trapped boundary within which we should carefully integrate fluxes
real_t BounceAverager::EvaluateBounceIntegral(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, function<real_t(real_t,real_t,real_t,real_t)> F, int_t *Flist){

    real_t xi0 = GetXi0(ir,i,j,fluxGridType);
    bool isTrapped = BounceSurfaceQuantity::IsTrapped(ir,i,j,fluxGridType,grid);
    real_t Bmin = fluxSurfaceAverager->GetBmin(ir,fluxGridType);
    real_t Bmax = fluxSurfaceAverager->GetBmax(ir,fluxGridType);

    real_t SingularPointCorrection = 1;

    // in inhomogeneous magnetic fields, for the cell centers, treat 
    // the trapped-passing boundary and origin more carefully
    // XXX: assumes P-Xi grid
    if(Bmin!=Bmax && fluxGridType==FLUXGRIDTYPE_DISTRIBUTION){
        real_t xi_f1 = GetXi0(ir,i,j,FLUXGRIDTYPE_P2);
        real_t xi_f2 = GetXi0(ir,i,j+1,FLUXGRIDTYPE_P2);
        real_t xiT = grid->GetRadialGrid()->GetXi0TrappedBoundary(ir);
        // if the cell contains, or is near, the trapped-passing boundaries, 
        // perform more careful (and computationally expensive) cell average
        real_t d1=fabs(xi_f1)-xiT;
        real_t d2=fabs(xi_f2)-xiT;
        bool containsTrappedBoundary = (xi_f1<xiT && xi_f2 >= xiT) || (xi_f1<=-xiT && xi_f2 >-xiT);
        bool isNearTrappedBoundary = min(fabs(d1),fabs(d2)) < NEAR_TRAPPED_BOUNDARY_THRESHOLD;
        if(containsTrappedBoundary || isNearTrappedBoundary){
            real_t p = grid->GetMomentumGrid(ir)->GetP1(i);
            return p*p*EvaluateCellAveragedBounceIntegralOverP2(ir,xi_f1,xi_f2,fluxGridType,F, Flist);
        
        // Treat singular xi0=0 case in inhomogeneous magnetic fields.
        // (unless the cell contains the trapping region, in which 
        //  case it would have been treated by the previous block)
        } else if(fabs(xi0)<100*realeps){
            // If cell center occurs at xi0=0, take the bounce integrals as 
            // the xi-average over the cell volume under the assumption that
            // the integrand varies linearly from its value at the upper cell face
            // to 0 at xi0=0.
            fluxGridType = FLUXGRIDTYPE_P2;
            j += 1;
            isTrapped = BounceSurfaceQuantity::IsTrapped(ir,i,j,fluxGridType,grid);
            xi0 = GetXi0(ir,i,j,fluxGridType);
            SingularPointCorrection = xi0/2;
        }
    }
    
    function<real_t(real_t,real_t,real_t,real_t)> F_eff;
    int_t *Flist_eff = nullptr;
    int_t Flist_copy[5];
    if(Flist != nullptr){
        for(len_t k=0;k<5;k++)
            Flist_copy[k] = Flist[k];
        Flist_eff = Flist_copy;
    }

    if (isTrapped){
        // trapped negative-pitch particles do not exist independently; their dynamics are described by the 
        // positive-pitch counterparts (since those are summed over both directions of motion). 
        if(xi0<0)
            return 0;
        // Sum quantity over both directions along the field line for trapped particle
        F_eff = [&](real_t x, real_t  y, real_t z, real_t w){return  (F(x,y,z,w) + F(-x,y,z,w)) ;};
        // i.e. if odd function in xi, set to zero, otherwise multiply by 2
        if(Flist != nullptr){
            if(Flist[0]%2==1)
                Flist_eff[0] = 0;
            else 
                Flist_eff[4] *= 2;
        }
    } else        
        F_eff = F;

    real_t theta_b1 = BounceSurfaceQuantity::Theta_B1(ir,i,j,fluxGridType,grid);
    real_t theta_b2 = BounceSurfaceQuantity::Theta_B2(ir,i,j,fluxGridType,grid);
    if(geometryIsSymmetric && theta_b1>0) // there are two identical and mirrored intervals which contribute
        SingularPointCorrection *= 2;
        
    real_t BounceIntegral = 0;

    bool integrateQAWS = false;
    // If using adaptive-integration setting, perform bounce integral with GSL quadrature
    if( ( (!isTrapped) && (integratePassingAdaptive)) || (isTrapped && integrateTrappedAdaptive) ){
        if(isTrapped && F_eff(0,1,1,1)!=0)
            integrateQAWS = true;
        gsl_function GSL_func; 
        real_t p;
        if(fluxGridType==FLUXGRIDTYPE_P1)
            p = grid->GetMomentumGrid(0)->GetP_f1(i,j);
        else if(fluxGridType==FLUXGRIDTYPE_P2)
            p = grid->GetMomentumGrid(0)->GetP_f2(i,j);
        else
            p = grid->GetMomentumGrid(0)->GetP(i,j);

        BounceIntegralParams params = {xi0, p, F, Flist_eff, ir, theta_b1, theta_b2, Bmin, fluxGridType, fluxSurfaceAverager, integrateQAWS}; 
        GSL_func.function = &(BounceIntegralFunction);
        GSL_func.params = &params;
        real_t epsabs = 0, epsrel = 1e-6, lim = gsl_adaptive->limit, error;
        if(integrateQAWS) // use QAWS if integrand is singular
            gsl_integration_qaws(&GSL_func, theta_b1, theta_b2, qaws_table,epsabs,epsrel,lim,gsl_adaptive,&BounceIntegral, &error);
        else
            gsl_integration_qag(&GSL_func, theta_b1, theta_b2,epsabs,epsrel,lim, QAG_KEY,gsl_adaptive,&BounceIntegral, &error);
        return SingularPointCorrection*BounceIntegral;
    }

    // otherwise continue and use the chosen fixed quadrature

    const real_t *BOverBmin = this->BOverBmin->GetData(ir,i,j,fluxGridType);
    const real_t *ROverR0   = this->ROverR0->GetData(ir,i,j,fluxGridType);
    const real_t *NablaR2   = this->NablaR2->GetData(ir,i,j,fluxGridType);
    const real_t *Metric    = this->Metric->GetData(ir,i,j,fluxGridType);

    len_t ntheta;
    const real_t *weights;
    real_t weightScaleFactor;
    if(isTrapped){
        ntheta  = ntheta_interp_trapped;
        weights = weights_trapped_ref;
        weightScaleFactor = theta_b2-theta_b1;
    } else {
        ntheta  = ntheta_interp_passing;
        weights = this->weights_passing;
        weightScaleFactor = 1;
    }

    real_t xiOverXi0,w;        
    for (len_t it = 0; it<ntheta; it++) {
        // treat the singular cylindrical case 
        if(BOverBmin[it]==1 || !xi0)
            xiOverXi0 = 1;
        else {
            real_t xi0Sq = xi0*xi0;
            xiOverXi0 = sqrt((1- BOverBmin[it] * (1-xi0Sq))/xi0Sq);
        }
        w = weightScaleFactor*weights[it];

        BounceIntegral += 2*M_PI*w*Metric[it]*F_eff(xiOverXi0,BOverBmin[it],ROverR0[it],NablaR2[it]);
    }        
    return SingularPointCorrection*BounceIntegral;    
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

    grid->SetBounceParameters(hasTrapped,
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
        real_t theta_Bmin = 0, theta_Bmax = 0;
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
                    FluxSurfaceAverager::FindBouncePoints(ir,Bmin, theta_Bmin, theta_Bmax, fluxSurfaceAverager, xi0, fluxGridType,
                                    &theta_b1[ir][pind],&theta_b2[ir][pind], gsl_fsolver, geometryIsSymmetric);
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
    BOverBmin->SetGridResolution(nr,np1,np2);
    ROverR0->SetGridResolution(nr,np1,np2);
    NablaR2->SetGridResolution(nr,np1,np2);
    Metric->SetGridResolution(nr,np1,np2);
}


/** 
 * Helper function for gsl integration that evaluates the bounce integral function
 */
struct PXiIntegralParams {len_t ir; fluxGridType fgType; function<real_t(real_t,real_t,real_t,real_t)> F; int_t *Flist; FluxSurfaceAverager *FSA;};
real_t BounceAverager::evalPXiIntegralFunc(real_t xi, void *par){
    PXiIntegralParams *params = (struct PXiIntegralParams*)par;
    len_t ir = params->ir;
    fluxGridType fluxGridType = params->fgType;
    function<real_t(real_t,real_t,real_t,real_t)> F = params->F;
    int_t *Flist = params->Flist;
    FluxSurfaceAverager *FSA = params->FSA;

    return FSA->EvaluatePXiBounceIntegralAtP(ir,xi,fluxGridType,F, Flist);
}


/**
 * Averages the bounce integral over xi from xi_l to xi_u 
 */
real_t BounceAverager::EvaluateCellAveragedBounceIntegralOverP2(len_t ir, real_t xi_l, real_t xi_u, fluxGridType fluxGridType, function<real_t(real_t,real_t,real_t,real_t)> F, int_t *F_list){
    real_t 
        partResult1 = 0,
        partResult2 = 0,
        partResult3 = 0;
    
    real_t dxi; 
    if(xi_u>xi_l)
        dxi = xi_u - xi_l;
    else if( fabs(xi_u-xi_l) < 100*realeps ) // simply evaluate the bounce integral at the point xi_u=xi_l
        return fluxSurfaceAverager->EvaluatePXiBounceIntegralAtP(ir,xi_u,fluxGridType,F,F_list);
    else 
        throw FVMException("BounceAverager: in EvaluateCellveragedBounceIntegral the upper xi value must be larger than, or equal to, the lower.");
    
    int key = GSL_INTEG_GAUSS41;
    real_t 
        epsabs = 0,
        epsrel = 1e-3,
        lim = gsl_adaptive->limit,
        error;
    
    real_t xiT = 0;
    if(fluxGridType == FLUXGRIDTYPE_RADIAL)
        xiT = grid->GetRadialGrid()->GetXi0TrappedBoundary_fr(ir);
    else 
        xiT = grid->GetRadialGrid()->GetXi0TrappedBoundary(ir);
    PXiIntegralParams params = {ir, fluxGridType, F, F_list, fluxSurfaceAverager}; 
    gsl_function gsl_func;
    gsl_func.function = &(evalPXiIntegralFunc);
    gsl_func.params = &params;
    
    // contribution from negative pitch passing region
    if(xi_l < -xiT){
        real_t pts[2] = {xi_l,-xiT};
        int npts = 2;
        if(xi_u < -xiT)
            gsl_integration_qag(&gsl_func,xi_l,xi_u,epsabs,epsrel,lim,key,gsl_adaptive,&partResult1,&error);
        else
            gsl_integration_qagp(&gsl_func,pts,npts,epsabs,epsrel,lim,gsl_adaptive,&partResult1,&error);
    }
    // contribution from positive pitch trapped region
    if(xi_u>0 && xi_l<xiT){
        if(xi_u<xiT){
            if(xi_l<0)
                gsl_integration_qag(&gsl_func,0,xi_u,epsabs,epsrel,lim,key,gsl_adaptive,&partResult2,&error);
            else 
                gsl_integration_qag(&gsl_func,xi_l,xi_u,epsabs,epsrel,lim,key,gsl_adaptive,&partResult2,&error);
        } else {
            if(xi_l<0){
                real_t pts[2] = {0,xiT};
                int npts = 2;
                gsl_integration_qagp(&gsl_func,pts,npts,epsabs,epsrel,lim,gsl_adaptive,&partResult2,&error);
            } else {
                real_t pts[2] = {xi_l,xiT};
                int npts = 2;
                gsl_integration_qagp(&gsl_func,pts,npts,epsabs,epsrel,lim,gsl_adaptive,&partResult2,&error);
            }
        }
    }
    // contribution from positive pitch passing region
    if(xi_u>xiT){
        real_t pts[2] = {xiT,xi_u};
        int npts = 2;
        if(xi_l<=xiT)
            gsl_integration_qagp(&gsl_func,pts,npts,epsabs,epsrel,lim,gsl_adaptive,&partResult3,&error);
        else 
            gsl_integration_qag(&gsl_func,xi_l,xi_u,epsabs,epsrel,lim,key,gsl_adaptive,&partResult3,&error);
    }    

    return (partResult1+partResult2+partResult3)/dxi;
}
