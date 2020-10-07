/**
 * Implementation of an overarching 'Grid' object.
 */

#include <algorithm>
#include <vector>
#include "FVM/Grid/Grid.hpp"


using namespace std;
using namespace DREAM::FVM;

/**
 * Constructor.
 */
Grid::Grid(RadialGrid *rg, MomentumGrid *mg, const real_t /*t0*/, FluxSurfaceAverager::quadrature_method qm_trapped, len_t ntheta_interp_trapped) {
    this->rgrid = rg;
    this->momentumGrids = new MomentumGrid*[rgrid->GetNr()];

    for (len_t i = 0; i < rgrid->GetNr(); i++)
        this->momentumGrids[i] = mg;


    FluxSurfaceAverager *FSA = rg->GetFluxSurfaceAverager();

    // if default ntheta, take same as in the flux surface averager 
    if(ntheta_interp_trapped==0)
        ntheta_interp_trapped = FSA->GetNTheta();

    bounceAverager = new BounceAverager(this, FSA,ntheta_interp_trapped,qm_trapped);
}

/**
 * Destructor.
 */
Grid::~Grid() {
    const len_t nr = this->GetNr();

    // Destroy momentum grids
    //   Since several, or even all, radii may share
    //   a single momentum grid, we should be careful
    //   not to try to double-free any momentum grid.
    vector<MomentumGrid*> deletedPtrs(nr);
    for (len_t i = 0; i < nr; i++) {
        MomentumGrid *p = this->momentumGrids[i];

        // Has the MomentumGrid been deleted already?
        if (find(deletedPtrs.begin(), deletedPtrs.end(), p) != deletedPtrs.end())
            continue;

        deletedPtrs.push_back(p);
        delete p;
    }
    
    DeallocateVprime();
    DeallocateBAvg();
    DeallocateBounceParameters();
    delete [] this->momentumGrids;
    delete this->rgrid;
    delete this->bounceAverager;

    if(avalancheDeltaHat != nullptr){
        for(len_t ir=0; ir<nr; ir++){
            delete [] avalancheDeltaHat[ir];
            delete [] avalancheDeltaHatNegativePitch[ir];
        }
        delete [] avalancheDeltaHat;
        delete [] avalancheDeltaHatNegativePitch;
    }
}

/*****************************
 * PUBLIC METHODS            *
 *****************************/
/**
 * Get the total number of cells on this grid,
 * including on the momentum grids at each radius.
 */
const len_t Grid::GetNCells() const {
    const len_t Nr = this->GetNr();
    len_t N = 0;

    for (len_t i = 0; i < Nr; i++)
        N += this->momentumGrids[i]->GetNCells();

    return N;
}

/**
 * Get the total number of cells on the radial flux grid,
 * including on the momentum grids at each radius.
 */
const len_t Grid::GetNCells_fr() const {
    const len_t Nr = this->GetNr();
    len_t N = 0;

    for (len_t i = 0; i < Nr; i++)
        N += this->momentumGrids[i]->GetNCells();

    // XXX here we assume that all momentum grids are the same
    // double count last momentum grid (since nr+1)
    N += this->momentumGrids[Nr-1]->GetNCells();

    return N;
}

/**
 * Get the total number of cells on the p1 flux grid.
 */
const len_t Grid::GetNCells_f1() const {
    const len_t Nr = this->GetNr();
    len_t N = 0;

    for (len_t i = 0; i < Nr; i++)
        N += this->momentumGrids[i]->GetNCells_f1();

    return N;
}

/**
 * Get the total number of cells on the p2 flux grid.
 */
const len_t Grid::GetNCells_f2() const {
    const len_t Nr = this->GetNr();
    len_t N = 0;

    for (len_t i = 0; i < Nr; i++)
        N += this->momentumGrids[i]->GetNCells_f2();

    return N;
}

/**
 * Integrate the given vector numerically over the entire
 * phase space (radius+momentum).
 *
 * vec: Vector to integrate numerically.
 */
real_t Grid::Integral(const real_t *vec) const {
    real_t I = 0;

    for (len_t ir = 0, offset=0; ir < this->GetNr(); ir++) {
        real_t VpVol = this->GetVpVol(ir);
        real_t dr    = this->GetRadialGrid()->GetDr(ir);
        I += this->IntegralMomentumAtRadius(ir, vec+offset) * VpVol * dr;

        offset += this->GetMomentumGrid(ir)->GetNCells();
    }

    return I;
}

/**
 * Integrate the given vector numerically over the momentum
 * grids at all radii. The result is stored in the vector
 * 'I', which must be of size 'nr'. The vector 'I' is
 * subsequently returned. If 'I' is a nullptr, a new vector
 * of size 'nr' is allocated and returned.
 *
 * vec: Vector to integrate numerically.
 * I:   Contains integral on return. Must be of size 'nr'. If
 *      'nullptr', this method allocates a new vector of size
 *      'nr', stores the result in it and returns it.
 */
real_t *Grid::IntegralMomentum(const real_t *vec, real_t *I) const {
    if (I == nullptr)
        I = new real_t[this->GetNr()];

    for (len_t ir = 0, offset=0; ir < this->GetNr(); ir++) {
        I[ir] = this->IntegralMomentumAtRadius(ir, vec+offset);
        offset += this->GetMomentumGrid(ir)->GetNCells();
    }

    return I;
}

/**
 * Integrate the given vector numerically over the momentum
 * grid corresponding to the specified radius.
 *
 * ir:  Index of radius to integrate over.
 * vec: Vector to integrate of size np1*np2.
 */
real_t Grid::IntegralMomentumAtRadius(const len_t ir, const real_t *vec) const {
    MomentumGrid *mg = this->GetMomentumGrid(ir);
    const len_t np1 = mg->GetNp1(), np2 = mg->GetNp2();
    const real_t *Vp = this->GetVp(ir);
    const real_t VpVol = this->GetVpVol(ir);

    real_t I = 0;
    for (len_t j = 0; j < np2; j++) {
        real_t dp2 = mg->GetDp2(j);

        for (len_t i = 0; i < np1; i++) {
            real_t dp1 = mg->GetDp1(i);
            len_t idx = j*np1 + i;

            I += vec[idx]*Vp[idx] * dp1*dp2;
        }
    }
    
    return I/VpVol;
}

/**
 * Rebuilds any non-static (i.e. time dependent) grids
 * used. This can be used if, for example, a dynamically
 * evolving magnetic equilibrium is used, or if some
 * grids are adaptive.
 *
 * t: Time to which re-build the grids for.
 */
bool Grid::Rebuild(const real_t t) {
    bool rgridUpdated = false, updated = false;

    // Re-build radial grid
    if (this->rgrid->NeedsRebuild(t))
        rgridUpdated = this->rgrid->Rebuild(t);

    updated = rgridUpdated;

    // Re-build momentum grids
    const len_t Nr = GetNr();
    for (len_t i = 0; i < Nr; i++) {
        if (this->momentumGrids[i]->NeedsRebuild(t, rgridUpdated))
            updated |= this->momentumGrids[i]->Rebuild(t, i, this->rgrid);
    }

    // Re-build jacobians
    if (updated)
        this->RebuildJacobians();

    return updated;
}

/**
 * Rebuilds magnetic-field data and initializes 
 * flux surface and bounce average calculations.
 */
void Grid::RebuildJacobians(){ 
    this->rgrid->RebuildJacobians(); 
    this->bounceAverager->Rebuild();
    RebuildBounceAveragedQuantities();
}

/**
 * Calculates and stores bounce averaged quantities.
 */
void Grid::RebuildBounceAveragedQuantities(){
 real_t 
    **BA_xi_f1,
    **BA_xi_f2, 
    **BA_xi2OverB_f1, 
    **BA_xi2OverB_f2,
    **BA_B3_f1,
    **BA_B3_f2,
    **BA_xi2B2_f1,
    **BA_xi2B2_f2,
    **BA_xiOverBR2;
    
    std::function<real_t(real_t,real_t,real_t,real_t)> F_xi = [](real_t xiOverXi0, real_t, real_t,real_t ){return xiOverXi0;};
    SetBounceAverage(BA_xi_f1, F_xi,FLUXGRIDTYPE_P1);
    SetBounceAverage(BA_xi_f2, F_xi,FLUXGRIDTYPE_P2);
    std::function<real_t(real_t,real_t,real_t,real_t)> F_xi2OverB = [](real_t xiOverXi0, real_t BOverBmin, real_t,real_t ){return xiOverXi0*xiOverXi0/BOverBmin;};
    SetBounceAverage(BA_xi2OverB_f1, F_xi2OverB,FLUXGRIDTYPE_P1);
    SetBounceAverage(BA_xi2OverB_f2, F_xi2OverB,FLUXGRIDTYPE_P2);
    std::function<real_t(real_t,real_t,real_t,real_t)> F_B3 = [](real_t , real_t BOverBmin, real_t,real_t ){return BOverBmin*BOverBmin*BOverBmin;};
    SetBounceAverage(BA_B3_f1, F_B3,FLUXGRIDTYPE_P1);
    SetBounceAverage(BA_B3_f2, F_B3,FLUXGRIDTYPE_P2);
    std::function<real_t(real_t,real_t,real_t,real_t)> F_xi2B2 = [](real_t xiOverXi0, real_t BOverBmin, real_t,real_t){return xiOverXi0*xiOverXi0*BOverBmin*BOverBmin;};
    SetBounceAverage(BA_xi2B2_f1, F_xi2B2,FLUXGRIDTYPE_P1);
    SetBounceAverage(BA_xi2B2_f2, F_xi2B2,FLUXGRIDTYPE_P2);
    std::function<real_t(real_t,real_t,real_t,real_t)> F_xiOverBR2 = [](real_t xiOverXi0, real_t BOverBmin, real_t ROverR0,real_t){return xiOverXi0/(BOverBmin*ROverR0*ROverR0);};
    SetBounceAverage(BA_xiOverBR2, F_xiOverBR2,FLUXGRIDTYPE_DISTRIBUTION);

    InitializeBAvg(BA_xi_f1,BA_xi_f2,BA_xi2OverB_f1, BA_xi2OverB_f2,BA_B3_f1,BA_B3_f2,
        BA_xi2B2_f1,BA_xi2B2_f2,BA_xiOverBR2);

    CalculateAvalancheDeltaHat();

}

/**
 * Calculate bounce average
 */
real_t Grid::CalculateBounceAverage(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F){
    return bounceAverager->CalculateBounceAverage(ir,i,j,fluxGridType,F);
}


/**
 * Calculate flux surface average
 */
real_t Grid::CalculateFluxSurfaceAverage(len_t ir, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t)> F){
    return rgrid->CalculateFluxSurfaceAverage(ir,fluxGridType,F);
}
    
/**
 * Helper method to set one bounce average
 */
void Grid::SetBounceAverage(real_t **&BA_quantity, std::function<real_t(real_t,real_t,real_t,real_t)> F, fluxGridType fluxGridType){
    len_t nr = GetNr() + (fluxGridType==FLUXGRIDTYPE_RADIAL);
    len_t np1, np2;
    BA_quantity = new real_t*[nr];
    for(len_t ir=0; ir<nr; ir++){
        MomentumGrid *mg = momentumGrids[ir];
        np1 = mg->GetNp1() + (fluxGridType==FLUXGRIDTYPE_P1);
        np2 = mg->GetNp2() + (fluxGridType==FLUXGRIDTYPE_P2);
        len_t ind_i0; // set to 1 if p(0,0)=0 since metric is singular
        if(fluxGridType==FLUXGRIDTYPE_P1)
            ind_i0 = (mg->GetP_f1(0,0)==0);
        else if(fluxGridType==FLUXGRIDTYPE_P2)
            ind_i0 = (mg->GetP_f2(0,0)==0);
        else 
            ind_i0 = (mg->GetP(0,0)==0);

        BA_quantity[ir] = new real_t[np1*np2];
        for(len_t j=0;j<np2;j++){
            if(ind_i0==1){
                real_t xi0;
                if(fluxGridType==FLUXGRIDTYPE_P1)
                    xi0 = mg->GetXi0_f1(0,j);
                else if(fluxGridType==FLUXGRIDTYPE_P1)
                    xi0 = mg->GetXi0_f1(0,j);
                else 
                    xi0 = mg->GetXi0(0,j);
                BA_quantity[ir][j*np1] = this->rgrid->CalculatePXiBounceAverageAtP(ir,0,xi0,fluxGridType,F);
            } 
            for(len_t i=ind_i0;i<np1;i++)
                BA_quantity[ir][j*np1+i] = CalculateBounceAverage(ir,i,j,fluxGridType,F);
        }
    }    
}

/**
 * Evaluates and stores the bounce- and cell averaged delta function 
 * that appears in the Rosenbluth-Putvinski avalanche source.
 * See documentation in doc/notes/theory for more details.
 */
void Grid::CalculateAvalancheDeltaHat(){
    if(avalancheDeltaHat != nullptr){
        for(len_t ir=0; ir<GetNr(); ir++){
            delete [] avalancheDeltaHat[ir];
            delete [] avalancheDeltaHatNegativePitch[ir];
        }
        delete [] avalancheDeltaHat;
        delete [] avalancheDeltaHatNegativePitch;
    }

    avalancheDeltaHat = new real_t*[GetNr()];
    avalancheDeltaHatNegativePitch = new real_t*[GetNr()];

    for(len_t ir=0; ir<GetNr(); ir++){
        MomentumGrid *mg = momentumGrids[ir];
        real_t VpVol = GetVpVol(ir);
        len_t np1 = GetNp1(ir);
        len_t np2 = GetNp2(ir);
        avalancheDeltaHat[ir] = new real_t[np1*np2];        
        avalancheDeltaHatNegativePitch[ir] = new real_t[np1*np2]; 
        for(len_t i=0; i<np1*np2; i++){
            avalancheDeltaHat[ir][i] = 0;
            avalancheDeltaHatNegativePitch[ir][i] = 0;
        }  
        for(len_t i=0; i<np1; i++)
            for(len_t j=0; j<np2; j++){
                real_t p = mg->GetP1(i);
                real_t xi_l = mg->GetP2_f(j);
                real_t xi_u = mg->GetP2_f(j+1);

                len_t j_tmp=j;
                // if negative-pitch trapped boundary, find index containing 
                // mirrored (-xi) cell to which we instead add the contribution
                if(IsTrapped_f2(ir,i,j+1) && xi_u<=0)
                    while(mg->GetP2_f(j_tmp+1)<-mg->GetP2(j) && j_tmp<np2)
                        j_tmp++;    

                // normalize contribution with dxi in new cell
                real_t fac = mg->GetDp2(j)/mg->GetDp2(j_tmp);
                real_t Vp = this->Vp[ir][j_tmp*np1+i];
                avalancheDeltaHat[ir][j_tmp*np1+i] += fac*bounceAverager->EvaluateAvalancheDeltaHat(ir,p,xi_l,xi_u,Vp, VpVol); \
                avalancheDeltaHatNegativePitch[ir][j_tmp*np1+i] += fac*bounceAverager->EvaluateAvalancheDeltaHat(ir,p,xi_l,xi_u,Vp, VpVol,-1); \
            }        
    }

}

/**
 * Set bounce averages
 */
void Grid::InitializeBAvg(
            real_t **xiAvg_f1, real_t **xiAvg_f2,
            real_t **xi2B2Avg_f1, real_t **xi2B2Avg_f2,
            real_t **B3_f1, real_t **B3_f2,
            real_t **xi2B2_f1, real_t **xi2B2_f2, real_t **xiOverBR2)
{
    DeallocateBAvg();
    this->BA_xi_f1                   = xiAvg_f1;
    this->BA_xi_f2                   = xiAvg_f2;
    this->BA_xi2OverB_f1             = xi2B2Avg_f1;
    this->BA_xi2OverB_f2             = xi2B2Avg_f2;
    this->BA_B3_f1                   = B3_f1;
    this->BA_B3_f2                   = B3_f2;
    this->BA_xi2B2_f1                = xi2B2_f1;
    this->BA_xi2B2_f2                = xi2B2_f2;
    this->BA_xiOverBR2               = xiOverBR2;   
}
/**
 * Deallocate bounce averages
 */
void Grid::DeallocateBAvg(){
    if (this->BA_xi_f1 == nullptr)
        return;
    
    for (len_t i = 0; i < GetNr(); i++) {
        delete [] this->BA_xi_f1[i];
        delete [] this->BA_xi_f2[i];
        delete [] this->BA_xi2OverB_f1[i];
        delete [] this->BA_xi2OverB_f2[i];
    }
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


/**
 * Set data for isTrapped and poloidal-angle bounce points 
 */
void Grid::SetBounceParameters(bool **isTrapped, bool **isTrapped_fr, 
            bool **isTrapped_f1, bool **isTrapped_f2, 
            real_t **theta_b1, real_t **theta_b1_fr, real_t **theta_b1_f1, real_t **theta_b1_f2, 
            real_t **theta_b2, real_t **theta_b2_fr, real_t **theta_b2_f1, real_t **theta_b2_f2 )
{
    DeallocateBounceParameters();
    this->isTrapped    = isTrapped;
    this->isTrapped_fr = isTrapped_fr;
    this->isTrapped_f1 = isTrapped_f1;
    this->isTrapped_f2 = isTrapped_f2;

    this->theta_b1    = theta_b1;
    this->theta_b1_fr = theta_b1_fr;
    this->theta_b1_f1 = theta_b1_f1;
    this->theta_b1_f2 = theta_b1_f2;

    this->theta_b2    = theta_b2;
    this->theta_b2_fr = theta_b2_fr;
    this->theta_b2_f1 = theta_b2_f1;
    this->theta_b2_f2 = theta_b2_f2;
}
/**
 * Deallocator
 */
void Grid::DeallocateBounceParameters(){
    if(isTrapped==nullptr)
        return;
    for(len_t ir=0; ir<GetNr(); ir++){
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
    for(len_t ir=0; ir<GetNr()+1; ir++){
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

/**
 * Set bounce-averaged metric Vprime
 */
void Grid::SetVp(real_t **Vp, real_t **Vp_fr, real_t **Vp_f1, real_t **Vp_f2, real_t **VpOverP2AtZero){
    DeallocateVprime();
    this->Vp = Vp;
    this->Vp_fr = Vp_fr;
    this->Vp_f1 = Vp_f1;
    this->Vp_f2 = Vp_f2;
    this->VpOverP2AtZero = VpOverP2AtZero;
}

/**
 * Deallocator
 */
void Grid::DeallocateVprime() {
    if (this->Vp == nullptr)
        return;

    for (len_t i = 0; i < GetNr(); i++) {
        delete [] this->Vp_f2[i];
        delete [] this->Vp_f1[i];
        delete [] this->Vp[i];
        delete [] this->VpOverP2AtZero[i];
    }
    for (len_t i = 0; i < this->GetNr()+1; i++)
        delete [] this->Vp_fr[i];

    delete [] this->Vp_f2;
    delete [] this->Vp_f1;
    delete [] this->Vp_fr;
    delete [] this->Vp;
    
    delete [] this->VpOverP2AtZero;
}
