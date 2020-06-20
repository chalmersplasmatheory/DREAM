/**
 * Implementation of the radial grid.
 */

#include <algorithm>
#include <vector>
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include "gsl/gsl_integration.h"


using namespace std;
using namespace DREAM::FVM;

/***********************
 * Constructors        *
 ***********************/
/**
 * Initialize an empty grid by only specifying the
 * grid size.
 *
 * rg: Object to use for (re-)generating the radial grid.
 * t0: Time to initialize grid at.
 * ntheta_interp: Poloidal angle resolution in quadrature for flux surface and bounce averages.
 */
RadialGrid::RadialGrid(RadialGridGenerator *rg, const real_t /*t0*/)
    : nr(rg->GetNr()), generator(rg) {

    // Build radial grid for the first time using
    // the given RadialGridGenerator
    // rg->Rebuild(t0, this);


}

/**
 * Destructor.
 */
RadialGrid::~RadialGrid() {
    //DeallocateMagneticField();
    
    DeallocateVprime();
    
    DeallocateGrid();
    DeallocateFSAvg();

    DeallocateMagneticData();

    delete this->generator;
}

/**
 * Deallocators.
 */
void RadialGrid::DeallocateGrid() {
    if (this->r == nullptr)
        return;

    delete [] this->dr_f;
    delete [] this->dr;
    delete [] this->r_f;
    delete [] this->r;
}

/*
void RadialGrid::DeallocateMagneticField() {
    if (this->theta_ref == nullptr)
        return;

    delete [] this->B_ref_f;
    delete [] this->B_ref;
    delete [] this->theta_ref;
    delete [] this->Bmin_ref;
    delete [] this->Bmin_ref_f;
    delete [] this->Bmax_ref;
    delete [] this->Bmax_ref_f;

}
*/


void RadialGrid::DeallocateVprime() {
    if (this->Vp == nullptr)
        return;

    for (len_t i = 0; i < this->nr; i++) {
        delete [] this->Vp_f2[i];
        delete [] this->Vp_f1[i];
        delete [] this->Vp[i];
        delete [] this->VpOverP2AtZero[i];
    }
    for (len_t i = 0; i < this->nr+1; i++)
        delete [] this->Vp_fr[i];

    delete [] this->Vp_f2;
    delete [] this->Vp_f1;
    delete [] this->Vp_fr;
    delete [] this->Vp;

    delete [] this->VpVol;
    delete [] this->VpVol_f;
    
    delete [] this->VpOverP2AtZero;
}



void RadialGrid::DeallocateFSAvg(){
    if (this->effectivePassingFraction == nullptr)
        return;
    
    for (len_t i = 0; i < this->nr; i++) {
        delete [] this->BA_xi_f1[i];
        delete [] this->BA_xi_f2[i];
        delete [] this->BA_xi2OverB_f1[i];
        delete [] this->BA_xi2OverB_f2[i];
    }

    delete [] this->FSA_B;
    delete [] this->FSA_B_f;
    delete [] this->FSA_B2;
    delete [] this->FSA_B2_f;
    delete [] this->effectivePassingFraction;
    delete [] this->FSA_nablaR2OverR2;
    delete [] this->FSA_nablaR2OverR2_f;
    delete [] this->FSA_1OverR2;
    delete [] this->FSA_1OverR2_f;
    delete [] this->BA_xi_f1;
    delete [] this->BA_xi_f2;
    delete [] this->BA_BOverBOverXi_f1;
    delete [] this->BA_BOverBOverXi_f2;
    delete [] this->BA_B3_f1;
    delete [] this->BA_B3_f2;
    delete [] this->BA_xi2B2_f1;
    delete [] this->BA_xi2B2_f2;
    delete [] this->BA_xiOverBR2;
}

/***************************
 * PUBLIC METHODS          *
 ***************************/
/**
 * Rebuilds any non-static (i.e. time dependent) grids
 * used. This can be used if, for example, a dynamically
 * evolving magnetic equilibrium is used, or if some
 * grids are adaptive.
 *
 * t: Time to which re-build the grids for.
 */
bool RadialGrid::Rebuild(const real_t t) {
    // Re-build radial grid
    if (this->generator->NeedsRebuild(t))
        return this->generator->Rebuild(t, this);
    else return false;
}

void RadialGrid::RebuildFluxSurfaceAveragedQuantities(MomentumGrid **momentumGrids){
 real_t 
    *effectivePassingFraction   = nullptr, 
    *effectivePassingFraction_f = nullptr, 
    *FSA_B2   = nullptr,
    *FSA_B2_f = nullptr,
    *FSA_B    = nullptr,
    *FSA_B_f  = nullptr,
    *FSA_nablaR2OverR2   = nullptr,
    *FSA_nablaR2OverR2_f = nullptr, 
    *FSA_1OverR2   = nullptr,
    *FSA_1OverR2_f = nullptr,
    **BA_xi_f1 = nullptr,
    **BA_xi_f2 = nullptr, 
    **BA_xi2OverB_f1 = nullptr, 
    **BA_xi2OverB_f2 = nullptr,
    **BA_B3_f1 = nullptr,
    **BA_B3_f2 = nullptr,
    **BA_xi2B2_f1 = nullptr,
    **BA_xi2B2_f2 = nullptr,
    **BA_xiOverBR2 = nullptr;
    
    std::function<real_t(real_t,real_t,real_t)> F_xi = [](real_t xiOverXi0, real_t, real_t ){return xiOverXi0;};
    SetBounceAverage(momentumGrids, BA_xi_f1, F_xi,FLUXGRIDTYPE_P1);
    SetBounceAverage(momentumGrids, BA_xi_f2, F_xi,FLUXGRIDTYPE_P2);
    std::function<real_t(real_t,real_t,real_t)> F_xi2OverB = [](real_t xiOverXi0, real_t BOverBmin, real_t ){return xiOverXi0*xiOverXi0/BOverBmin;};
    SetBounceAverage(momentumGrids, BA_xi2OverB_f1, F_xi2OverB,FLUXGRIDTYPE_P1);
    SetBounceAverage(momentumGrids, BA_xi2OverB_f2, F_xi2OverB,FLUXGRIDTYPE_P2);
    std::function<real_t(real_t,real_t,real_t)> F_B3 = [](real_t , real_t BOverBmin, real_t ){return BOverBmin*BOverBmin*BOverBmin;};
    SetBounceAverage(momentumGrids, BA_B3_f1, F_B3,FLUXGRIDTYPE_P1);
    SetBounceAverage(momentumGrids, BA_B3_f2, F_B3,FLUXGRIDTYPE_P2);
    std::function<real_t(real_t,real_t,real_t)> F_xi2B2 = [](real_t xiOverXi0, real_t BOverBmin, real_t){return xiOverXi0*xiOverXi0*BOverBmin*BOverBmin;};
    SetBounceAverage(momentumGrids, BA_xi2B2_f1, F_xi2B2,FLUXGRIDTYPE_P1);
    SetBounceAverage(momentumGrids, BA_xi2B2_f2, F_xi2B2,FLUXGRIDTYPE_P2);
    std::function<real_t(real_t,real_t,real_t)> F_xiOverBR2 = [](real_t xiOverXi0, real_t BOverBmin, real_t ROverR0){return xiOverXi0/(BOverBmin*ROverR0*ROverR0);};
    SetBounceAverage(momentumGrids, BA_xiOverBR2, F_xiOverBR2,FLUXGRIDTYPE_DISTRIBUTION);
    
    
    SetFluxSurfaceAverage(FSA_1OverR2,FSA_1OverR2_f, [](real_t , real_t ROverR0, real_t ){return 1/(ROverR0*ROverR0);} );
    SetFluxSurfaceAverage(FSA_B,FSA_B_f, [](real_t BOverBmin, real_t , real_t ){return BOverBmin;} );
    SetFluxSurfaceAverage(FSA_B2,FSA_B2_f, [](real_t BOverBmin, real_t , real_t ){return BOverBmin*BOverBmin;} );
    SetFluxSurfaceAverage(FSA_nablaR2OverR2,FSA_nablaR2OverR2_f, [](real_t , real_t ROverR0, real_t NablaR2){return NablaR2/(ROverR0*ROverR0);} );
    

    SetEffectivePassingFraction(effectivePassingFraction,effectivePassingFraction_f, FSA_B2, FSA_B2_f);


    InitializeFSAvg(effectivePassingFraction,effectivePassingFraction_f,
        FSA_B,FSA_B_f,FSA_B2,FSA_B2_f,FSA_1OverR2, FSA_1OverR2_f,FSA_nablaR2OverR2,FSA_nablaR2OverR2_f, 
        BA_xi_f1,BA_xi_f2,BA_xi2OverB_f1, BA_xi2OverB_f2,BA_B3_f1,BA_B3_f2,
        BA_xi2B2_f1,BA_xi2B2_f2,BA_xiOverBR2);

}


void RadialGrid::SetFluxSurfaceAverage(real_t *&FSA_quantity, real_t *&FSA_quantity_f, std::function<real_t(real_t,real_t,real_t)> F){
    FSA_quantity   = new real_t[GetNr()];
    FSA_quantity_f = new real_t[GetNr()+1];

    for(len_t ir=0; ir<GetNr(); ir++){
        FSA_quantity[ir] = CalculateFluxSurfaceAverage(ir, FLUXGRIDTYPE_DISTRIBUTION, F);
    }

    for(len_t ir=0; ir<GetNr()+1; ir++){
        FSA_quantity_f[ir] = CalculateFluxSurfaceAverage(ir, FLUXGRIDTYPE_RADIAL, F);
    }
}
void RadialGrid::SetBounceAverage(MomentumGrid **momentumGrids, real_t **&BA_quantity, std::function<real_t(real_t,real_t,real_t)> F, fluxGridType fluxGridType){
    len_t nr = GetNr() + (fluxGridType==FLUXGRIDTYPE_RADIAL);
    len_t np1, np2;
    BA_quantity = new real_t*[nr];
    for(len_t ir=0; ir<nr; ir++){
        MomentumGrid *mg = momentumGrids[ir];
        np1 = mg->GetNp1() + (fluxGridType==FLUXGRIDTYPE_P1);
        np2 = mg->GetNp2() + (fluxGridType==FLUXGRIDTYPE_P2);
        len_t ind_i0; // set to 1 if p(0,0)=0 since bounce average is singular
        if(fluxGridType==FLUXGRIDTYPE_P1)
            ind_i0 = (mg->GetP_f1(0,0)==0);
        else if(fluxGridType==FLUXGRIDTYPE_P1)
            ind_i0 = (mg->GetP_f2(0,0)==0);
        else 
            ind_i0 = (mg->GetP(0,0)==0);

        BA_quantity[ir] = new real_t[np1*np2];
        for(len_t i=ind_i0;i<np1;i++)
            for(len_t j=0;j<np2;j++)
                BA_quantity[ir][j*np1+i] = CalculateBounceAverage(mg,ir,i,j,fluxGridType,F);
    } 
}

/*
void RadialGrid::SetBounceAverage(MomentumGrid **momentumGrids, real_t **&BA_quantity_f1, real_t **&BA_quantity_f2, std::function<real_t(real_t,real_t)> F){
    BA_quantity_f1 = new real_t*[GetNr()];
    BA_quantity_f2 = new real_t*[GetNr()];
    FVM::fluxGridType fluxGridType;

    for(len_t ir=0; ir<GetNr(); ir++){
        MomentumGrid *mg = momentumGrids[ir];
        len_t np1 = mg->GetNp1();
        len_t np2 = mg->GetNp2();
        BA_quantity_f1[ir] = new real_t[(np1+1)*np2];
        BA_quantity_f2[ir] = new real_t[np1*(np2+1)];

        fluxGridType = FVM::FLUXGRIDTYPE_P1;
        // TODO: Now ignore setting bounce averages at p=0. Could fix in the future, but need to revise GetMetric (and return sqrt(g)/p^2 instead)
        for (len_t i = 1; i<np1+1; i++){
            for (len_t j=0; j<np2; j++){
                BA_quantity_f1[ir][j*(np1+1)+i] = CalculateBounceAverage(mg,  ir,  i,  j,  fluxGridType, F);
            }
        }

        fluxGridType = FVM::FLUXGRIDTYPE_P2;
        for (len_t i = 0; i<np1; i++){
            for (len_t j=0; j<np2+1; j++){
                BA_quantity_f2[ir][j*np1+i] = CalculateBounceAverage(mg,  ir,  i,  j,  fluxGridType, F);
            }
        }  
    }     
    
}
*/


/**
 * The function is used in the evaluation of the effective passing fraction,
 * and represents x / <1-x B/Bmax>
 */
struct EPF_params {real_t BminOverBmax; len_t ir; RadialGrid *rGrid; fluxGridType fgType;};
real_t RadialGrid::effectivePassingFractionIntegrand(real_t x, void *p){
    struct EPF_params *params = (struct EPF_params *) p;
    RadialGrid *rGrid = params->rGrid;
    real_t BminOverBmax = params->BminOverBmax; 
    len_t ir = params->ir;
    fluxGridType fluxGridType = params->fgType;
    std::function<real_t(real_t,real_t,real_t)> fluxAvgFunc = [x,BminOverBmax](real_t BOverBmin,real_t, real_t){
        return sqrt(1 - x * BminOverBmax * BOverBmin );
    };
    return x/ rGrid->CalculateFluxSurfaceAverage(ir, fluxGridType, fluxAvgFunc);
}

void RadialGrid::SetEffectivePassingFraction(real_t *&EPF, real_t *&, real_t *FSA_B2, real_t*){
    gsl_integration_workspace *gsl_w = gsl_integration_workspace_alloc(1000);
    gsl_function EPF_func;
    EPF_params paramstruct;
    real_t EPF_integral;
    EPF = new real_t[GetNr()];
    real_t error = 0;
    for (len_t ir=0; ir<GetNr(); ir++){
        real_t Bmin = GetBmin(ir);
        real_t Bmax = GetBmax(ir);
        real_t BminOverBmax = Bmin/Bmax;
        paramstruct = {BminOverBmax,ir,this,FLUXGRIDTYPE_DISTRIBUTION};
        EPF_func.function = &(effectivePassingFractionIntegrand);
        EPF_func.params = &paramstruct;
        real_t epsabs = 0, epsrel = 1e-4, lim = gsl_w->limit;
        gsl_integration_qags(&EPF_func, 0,1,epsabs,epsrel,lim,gsl_w,&EPF_integral, &error);
        EPF[ir] = (3.0/4) * Bmin*Bmin*FSA_B2[ir] / (Bmax*Bmax) * EPF_integral;
    }

    // We will probably not need this most of the time, so commenting it out for now
    /* 
    bool rFluxGrid = true;
    for (len_t ir=0; ir<GetNr()+1; ir++){
        real_t Bmin = GetBmin_f(ir);
        real_t Bmax = GetBmax_f(ir);
        real_t BminOverBmax = Bmin/Bmax;
        EPF_params paramstruct = {BminOverBmax,ir,mgnQtyHandler,rFluxGrid}; 
        EPF_func.function = &(effectivePassingFractionIntegrand);
        EPF_func.params = &paramstruct;
        gsl_integration_qags(&EPF_func, 0,1,0,1e-7,1000,gsl_w,&EPF_integral,nullptr );
        EPF_f[ir] = (3/4) * GetFSA_B2_f(ir) / (Bmax*Bmax) * EPF_integral;
    }
    */

    gsl_integration_workspace_free(gsl_w);
    
}

/**
 * Sets reference magnetic field quantities that have been
 * generated by a RadialGridGenerator (FluxSurfaceAverager
 * takes ownership of these and will deallocate them).
 */
void RadialGrid::SetReferenceMagneticFieldData(
    len_t ntheta_ref, real_t *theta_ref,
    real_t **B_ref, real_t **B_ref_f,
    real_t **Jacobian_ref, real_t **Jacobian_ref_f,
    real_t **ROverR0_ref ,real_t **ROverR0_ref_f, 
    real_t **NablaR2_ref, real_t **NablaR2_ref_f,
    real_t *Bmin, real_t *Bmin_f,
    real_t *Bmax, real_t *Bmax_f,
    real_t *BtorGOverR0, real_t *BtorGOverR0_f,
    real_t R0
){
    DeallocateMagneticData();

    this->Bmin           = Bmin;
    this->Bmin_f         = Bmin_f;
    this->Bmax           = Bmax;
    this->Bmax_f         = Bmax_f;
    this->BtorGOverR0    = BtorGOverR0;
    this->BtorGOverR0_f  = BtorGOverR0_f;
    this->R0             = R0;

    fluxSurfaceAverager->SetReferenceMagneticFieldData(
        ntheta_ref, theta_ref, B_ref, B_ref_f,
        Jacobian_ref, Jacobian_ref_f,
        ROverR0_ref , ROverR0_ref_f, 
        NablaR2_ref,  NablaR2_ref_f
    );
}


real_t RadialGrid::CalculateFluxSurfaceAverage(len_t ir, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t)> F){
    return fluxSurfaceAverager->CalculateFluxSurfaceAverage(ir, fluxGridType, F);
}

