/**
 * Implementation of a diffusive transport term which can be applied to both
 * kinetic and fluid grids, and which allows one to prescribe the diffusion
 * coefficient in time and phase space.
 */

#include <type_traits>
#include "DREAM/Equations/Fluid/SvenssonTransport.hpp"

/**
 * Constructor.
 */
template<typename T>
DREAM::SvenssonTransport<T>::SvenssonTransport(
    DREAM::FVM::Grid *grid,
    real_t pStar,
    DREAM::FVM::UnknownQuantityHandler *unknowns,
    DREAM::RunawayFluid *REFluid,
    FVM::Interpolator3D *interp3d,
    enum FVM::Interpolator3D::momentumgrid_type mtype
) : T(grid),
    nr(grid->GetNr()), np(grid->GetNp1(0)),// np(interp3d->GetNx3()),
    EID(unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD)),
    pStar(pStar),
    unknowns(unknowns), REFluid(REFluid),
    interp3d(interp3d), mtype(mtype)
{
    //this->EID = this->unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD); 

    printf("pStar = %f\n",pStar); // DEBUG
    
    /**
     * YYY
     * Do we need a momentum grid, or can we not just use the values
     * provided to interp3d, for the integration (without doing any
     * interpolation?
     */
    
    this->integrand = new real_t[this->np];
    
    //this->p = interp3d->GetX3();
    this->p = this->grid->GetMomentumGrid(0)->GetP1();

    // GSL integral (FVM/Grid/BounceAverager)
    // PA average coeffs...
    // (Begin at only xi=1)
    // Do the averaging to the x3
    this->coeff = new real_t[(nr+1)*np];


    len_t nxi=this->grid->GetNp2(0);
    // DEBUG
    // printf("nr=%lu,\tnp=%lu,\tnxi=%lu\n\n",nr,np,nxi);
    // fflush(stdout);

    real_t *out = new real_t[(nr+1)*np*nxi];
    interp3d->Eval(this->grid, mtype, FVM::FLUXGRIDTYPE_RADIAL, out);

    // DEBUG:
    std::cout << "Current momentumgrid_type: " << mtype ;
    printf("\n");
    std::cout << "Reference grid type (pxi): "<< FVM::Interpolator3D::momentumgrid_type::GRID_PXI;
    printf("\n");
    

    real_t avg, xiRange;
    const real_t *dxi = this->grid->GetMomentumGrid(0)->GetDp2();
    for (len_t ir=0, offset=0; ir < nr+1 ; ir++){
        printf("ir = %2lu",ir); // DEBUG
        for (len_t i=0; i < this->np ; i++){
            // Do the GSL integration for PA averaging
            avg=0;
            xiRange=0;
            for (len_t j=0; j<nxi; j++){
                //avg+= out[ir*nr+offset+j] * dxi[j];
                // printf("%f, ",out[(ir*nxi+j)*np + i]); // DEBUG
                avg += out[(ir*nxi+j)*np + i]*dxi[j];
                xiRange += dxi[j];
                // YYY Note that the distribution function is still missing!
            }
            coeff[i+offset] = avg/xiRange;//GSL stuff;
            printf("\t\tip = %2lu, coeff=%f, ", i, avg); // DEBUG
        }
        offset+=this->np;
        printf("\n"); fflush(stdout); // DEBUG
    }
    
    printf("\n"); fflush(stdout); // DEBUG
    
    delete [] out;
}

/**
 * Destructor.
 */
template<typename T>
DREAM::SvenssonTransport<T>::~SvenssonTransport() {
    delete [] this->integrand;
}



/**
 * Rebuild this term by evaluating and setting the diffusion
 * coefficient for the next time step.
 */
template<typename T>
void DREAM::SvenssonTransport<T>::Rebuild(
    const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler*
    ) {
    //const real_t *c = this->prescribedCoeff->Eval(t);
    
    //const len_t nr = this->grid->GetNr();
    
    //real_t *dp;
    
    
    // Iterate over the radial flux grid...
    for (len_t ir = 0; ir < this->nr+1; ir++) {
        
        // The varaible to be added to
        const real_t *dp = this->grid->GetMomentumGrid(0)->GetDp1();
        real_t pIntCoeff = 0;
        
        //const real_t *integrandArray = this->EvaluateIntegrand(ir);
        this->EvaluateIntegrand(ir);
        
        for (len_t i = 0; i < this->np; i++) {
            // The actual integration in p
            pIntCoeff += this->integrand[i] * dp[i];
                // YYY Jacobian??? * this->grid->GetVp(ir,i,0); 
        }

        this->_setcoeff(ir, pIntCoeff);        
        //offset += this->np;
    }
}



/**
 * Helper function for calculating the inverse of p-bar, with the
 * (optional) additional calculation of the derivative of
 * p-bar-inverse. 
 *
 * `p-bar` is the name given to the factor dividing `-(p - p*)` in the
 * exponential of eqn (4.2) in Svensson et al. 2020
 * [https://arxiv.org/abs/2010.07156v1].
 * 
 * These values are calculated on the flux grid, meaning that
 * interpolation (and extrapolation) from the cell grid is being
 * performed. This is done via inter-/extrapolation of p-bar-inverse,
 * instead of first inter-/extrsapolating the values going into p-bar.
 */
template<typename T>
real_t DREAM::SvenssonTransport<T>::GetPBarInv_f(len_t ir, real_t *dr_pBarInv_f){
    // Need interpolation from cell grid to flux grid:
    // pBar_f[0]=pBar[0]
    // pBar_f[ir]  = (pBar[ir-1] + pBar[ir] )*.5
    // pBar_f[nr]= extrapolate


    // Inverse of p-bar on the Flux grid, with additional helper variable.
    real_t pBarInv_f, tmp_pBarInv_f; 

    // Essential values taken on the raidal (cell) grid
    const real_t *E = this->unknowns->GetUnknownData(this->EID);
    const real_t *EcEff = this->REFluid->GetEffectiveCriticalField();
    const real_t *tauRel = this->REFluid->GetElectronCollisionTimeRelativistic();
    const real_t *gamma_r = this->REFluid->GetAvalancheGrowthRate();

    // Grid step size in the radial grid for the derivative.
    const real_t *dr_f = this->grid->GetRadialGrid()->GetDr_f();
    const real_t *dr = this->grid->GetRadialGrid()->GetDr(); 

    
    // Interpolating (extrapolating) the inverse of p bar onto the
    // flux grid.
    if (ir == 0) {
        // Zero flux at r = 0. Therefore choose the value at "ir=1/2".
        pBarInv_f = tauRel[0] * gamma_r[0] / (E[0]-EcEff[0]);
        
        if(dr_pBarInv_f != nullptr)
            *dr_pBarInv_f = 0.0;
    }
    else if (ir == this->nr) {
        // Linearly extrapolating the value at the end point from the
        // two previous points.
        //
        // N.B.! The extrapolation assume that the grid cell size is
        // uniform, and that the extrapolated value lies half a grid
        // cell away from the last point.

        // pBarInv_f  = 1.5 * (tauRel[ir-1] * gamma_r[ir-1] / (E[ir-1]-EcEff[ir-1]));
        // pBarInv_f -= 0.5 * (tauRel[ir-2] * gamma_r[ir-2] / (E[ir-2]-EcEff[ir-2]));
        pBarInv_f     = tauRel[ir-1] * gamma_r[ir-1] / (E[ir-1]-EcEff[ir-1]);
        tmp_pBarInv_f = tauRel[ir-2] * gamma_r[ir-2] / (E[ir-2]-EcEff[ir-2]);

        // N.B.! This order of operations is important
        
        if(dr_pBarInv_f != nullptr){
            // Derivative:
            *dr_pBarInv_f = (pBarInv_f - tmp_pBarInv_f) / dr_f[ir-2];
            // Extrapolation:
            pBarInv_f += (*dr_pBarInv_f) * 0.5*dr[ir-1];
        }
        else{
        pBarInv_f += (pBarInv_f - tmp_pBarInv_f) * 0.5*dr[ir-1]/dr_f[ir-2];
        }
        // The above is the same as:
        // pBarInv_f = tmp_pBarInv_f
        //     + (pBarInv_f - tmp_pBarInv_f)
        //     * (0.5*dr[ir-2] + dr[ir-1])/dr_f[ir-2];


        // This is for uniform step sizes:
        // pBarInv_f *= 1.5;
        // pBarInv_f -= 0.5 * tmp_pBarInv_f; 
    }
    else {
        // In the middle, we simply linearly interpolate
        tmp_pBarInv_f = tauRel[ir-1] * gamma_r[ir-1] / (E[ir-1]-EcEff[ir-1]);
        pBarInv_f  = tauRel[ir] * gamma_r[ir] / (E[ir]-EcEff[ir]);

        // N.B.! This order of operations is important!
        if(dr_pBarInv_f != nullptr){
            // Derivative:
            *dr_pBarInv_f = (pBarInv_f - tmp_pBarInv_f) / dr_f[ir-1]; 
            // Interpolation:
            //pBarInv_f = tmp_pBarInv_f + (*dr_pBarInv_f) * 0.5*dr[ir-1];
            pBarInv_f -= (*dr_pBarInv_f) * 0.5*dr[ir];
        }
        else{
        pBarInv_f = ( dr[ir]*tmp_pBarInv_f + dr[ir-1]*pBarInv_f ) * 0.5 / dr_f[ir-1];
        }
        // The above is the same as: 
        // pBarInv_f = tmp_pBarInv_f + (pBarInv_f - tmp_pBarInv_f) * 0.5*dr[ir-1]/dr_f[ir-1]

        // This is for uniform step sizes:
        // pBarInv_f += tmp_pBarInv_f;
        // pBarInv_f *= 0.5;
    }

    return pBarInv_f;
}
