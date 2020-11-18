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
    struct dream_4d_data *inputData4dStruct,
    enum FVM::Interpolator1D::interp_method timeInterpMethod 
) : T(grid),
    nr_f(grid->GetNr()+1),
    nt(inputData4dStruct->nt), nr(inputData4dStruct->nr),
    np1(inputData4dStruct->np1), np2(inputData4dStruct->np2),
    EID(unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD)),
    pStar(pStar),
    t(inputData4dStruct->t), r(inputData4dStruct->r),
    p1(inputData4dStruct->p1), p2(inputData4dStruct->p2),
    coeff4dInput(inputData4dStruct->x), // Size nt-by-(nr*np2*np1)
    inputMomentumGridType(inputData4dStruct->gridtype),
    inputInterp3dMethod(inputData4dStruct->ps_interp),
    unknowns(unknowns), REFluid(REFluid),
    timeInterpMethod(timeInterpMethod)
{
    // Checks that pStar is valid.
    if ( pStar <= 0 ) {
        throw DREAMException(
            "Invalid lower momentum bound, pStar=%.2f, must be >0.",
            pStar
            );
    }
    // Checks that the input data is defined on a p-xi momentum grid.
    if (inputMomentumGridType != FVM::Interpolator3D::GRID_PXI ){
        // SvenssonTransport is, at the moment, only defined for p-xi.
        throw DREAMException(
            "Wrong input momentum grid, input must be specified on p-xi grid."
            );
    }
    else{
        xi=p2;
        nxi=np2;
        // YYY the lower bound of pstar has not yet been implemented
        // YYY set p from pstar and up, p[0]=pstar        
        p=p1;

        np=np1; // YYY change later
    }


    this->coeffTRXiP = new real_t[ nt * nr_f * np * nxi ];
    this->coeffRP = new real_t[ nr_f * np ];
    this->integrand = new real_t[ np ];
    
    InterpolateCoefficient();
}


/**
 * Destructor.
 */
template<typename T>
DREAM::SvenssonTransport<T>::~SvenssonTransport() {
    // YYY Delete input data!
    delete [] this->coeffTRXiP;
    delete [] this->coeffRP;
    delete [] this->integrand;
    if (this->interpTCoeff != nullptr) {
        delete this->interpTCoeff;
    }
}




/**
 * Interpolate the time-depentent input coefficient onto the DREAM
 * radial flux-grid, then storing that data in `interpTCoeff`, which
 * is then used to evaluate the time dependance of the coefficients.
 *
 * This function creates a 1D interpolator, which is then used to
 * evaluate the time dependence of the coefficients.
 */
template<typename T>
void DREAM::SvenssonTransport<T>::InterpolateCoefficient() {    
    const len_t N  =  this->nr_f * this->nxi * this->np;
    
    for (len_t it = 0, offset = 0; it < nt; it++) {
        // YYY How do we ensure that we actually get p-xi and not ppar-pperp?
        DREAM::FVM::Interpolator3D intp3d_tmp(
            nr, np2, np1, r, p2, p1, coeff4dInput[it],
            inputMomentumGridType, inputInterp3dMethod, false
            );
        // Interpolating the coefficients onto the r_f grid used by
        // DREAM. We also make sure that the data is converted to a
        // xi-p grid.
        intp3d_tmp.Eval(nr_f, np2, np1, this->grid->GetRadialGrid()->GetR_f(), xi,p,
                        FVM::Interpolator3D::momentumgrid_type::GRID_PXI, coeffTRXiP+offset);
        offset+=N;
    }
    // Note that `coeffTRXiP` now contains r_f, xi and p data for _every_ time step.

    if (this->interpTCoeff != nullptr) {
        delete this->interpTCoeff;
    }
    // YYY Put the input 4d data into here and do the interpolation onto r_f in xiAverage!
    this->interpTCoeff = new DREAM::FVM::Interpolator1D(nt, N, t, coeffTRXiP, timeInterpMethod);
}





template<typename T>
void DREAM::SvenssonTransport<T>::xiAverage(const real_t *c){
    // Input `c` is the r-xi-p coefficient data (nr_f * nxi * np)
    // Writing xi-averaged data to `this->coeffRP` (nr_f * np)

    
    printf("         |");                                // DEBUG
    printf("           |");                              // DEBUG
    for(len_t i=0; i<np;i++) printf("  ip |  coeff |"); // DEBUG
    printf("\n");                                        // DEBUG
    
    if (nxi > 1) { 
        for (len_t ir=0, offset=0; ir < nr_f ; ir++){ // for radius
            printf("ir = %3lu | ",ir); // DEBUG
            printf("r = %0.3f | ",this->grid->GetRadialGrid()->GetR_f()[ir]); // DEBUG
    
            for (len_t i=0; i < this->np ; i++){ // for momentum
                // Do the GSL integration for PA averaging
                real_t avg = 0; // Varaible containing the xi average
                // GSL integral (example in `FVM/Grid/BounceAverager`)
                for (len_t j=0, j_offset=0; j<nxi-1; j++){ // for xi
                    len_t ind = offset*nxi + j_offset + i;
                    // avg += 0.5 * (c[ind] + c[ind+np]) * (xi[j+1]-xi[j]);
                    
                    real_t p_sq = this->p[i] * this->p[i];
                    real_t E = this->EvalOnFluxGrid(ir,
                        this->unknowns->GetUnknownData(this->EID));
                    real_t Zeff = this->EvalOnFluxGrid(ir,
                        this->REFluid->GetIonHandler()->evaluateZeff());
                    real_t w = 2.0 * E * p_sq
                        / ( (1.+Zeff) * sqrt(1+p_sq) );
                    // printf("w = %f\n",w);fflush(stdout); // DEBUG
                    real_t f1=c[ind] * exp(-w * this->xi[j]);
                    real_t f2=c[ind+np] * exp(-w * this->xi[j+1]);
                    avg += 0.25 * w / sinh(w) * (f1 + f2) * (xi[j+1]-xi[j]);
                    // The prefactor is 0.25 since there is one 0.5
                    // from the actual integral, and 0.5 from the
                    // trapz method.

                    j_offset += np;
                }
                this->coeffRP[i+offset] = avg / (xi[nxi-1]-xi[0]);
                printf("%3lu | %0.4f | ", i, this->coeffRP[i+offset]); // DEBUG
            }
            offset+=this->np;
            printf("\n"); fflush(stdout); // DEBUG
        }
    }
    else{
        for (len_t ir=0, offset=0; ir < nr_f ; ir++){ // for radius
            for (len_t i=0; i < this->np ; i++){ // for momentum
                this->coeffRP[i+offset] = c[offset + i];
            }
            offset+=this->np;
        }
    }
    
    // printf("%f\n",coeffRXiP[nr_f*nxi*np-1]); // DEBUG
    // printf("\n"); fflush(stdout);   // DEBUG
    
    // delete [] coeffRXiP;    
}



/**
 * Rebuild this term by evaluating and setting the advection and
 * diffusion terms for the next time step.
 */
template<typename T>
void DREAM::SvenssonTransport<T>::Rebuild(
    const real_t t, const real_t, DREAM::FVM::UnknownQuantityHandler*
    ) {

    // printf("t = %0.2f\n",t); fflush(stdout); // DEBUG
    
    const real_t *c = this->interpTCoeff->Eval(t);

    // Note that the average has to be redone for every timestep,
    // since the prescribed distributionfunction changes with time.
    xiAverage(c);
    
    if (np > 1){
        // Iterate over the radial flux grid...
        for (len_t ir = 0; ir < this->nr_f; ir++) {
            
            // The varaible to be added to
            // const real_t *dp = this->grid->GetMomentumGrid(0)->GetDp1();
            real_t pIntCoeff = 0;
            
            this->EvaluateIntegrand(ir);
            
            // The actual integration in p
            for (len_t i = 0; i < this->np-1; i++) {
                //pIntCoeff += this->integrand[i] * dp[i];
                pIntCoeff += 0.5 * (this->integrand[i] + this->integrand[i+1])
                    * (this->p[i+1] - this->p[i]);
                // YYY Jacobian??? * this->grid->GetVp(ir,i,0); 
            }
            //printf("pIntCoeff = %f\n",pIntCoeff);  fflush(stdout); // DEBUG
            this->_setcoeff(ir, pIntCoeff);
        }
    }
    else{
        for (len_t ir = 0; ir < this->nr_f; ir++) {
            // `coeffRP` is of size nr_f*np
            // printf("pIntCoeff = %f\n",coeffRP[ir]);  fflush(stdout); // DEBUG
            this->_setcoeff(ir, coeffRP[ir]);
        }
    }
}











/**
 * YYY Consider using EvalOnFluxGrid instead of almost the same
 * function for pBarInv.
 * 
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
    const real_t *E       = this->unknowns->GetUnknownData(this->EID);
    const real_t *EcEff   = this->REFluid->GetEffectiveCriticalField();
    const real_t *tauRel  = this->REFluid->GetElectronCollisionTimeRelativistic();
    const real_t *gamma_r = this->REFluid->GetAvalancheGrowthRate();

    // Grid step size in the radial grid for the derivative.
    const real_t *dr_f    = this->grid->GetRadialGrid()->GetDr_f();
    const real_t *dr      = this->grid->GetRadialGrid()->GetDr(); 

    
    // Interpolating (extrapolating) the inverse of p bar onto the
    // flux grid.
    if (ir == 0) {
        // Zero flux at r = 0. Therefore choose the value at "ir=1/2".
        pBarInv_f = tauRel[0] * gamma_r[0] / (E[0]-EcEff[0]);
        
        if(dr_pBarInv_f != nullptr)
            *dr_pBarInv_f = 0.0;
    }
    else if (ir == this->nr_f - 1) {
        // Linearly extrapolating the value at the end point from the
        // two previous points.

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
            //pBarInv_f = ( dr[ir]*tmp_pBarInv_f + dr[ir-1]*pBarInv_f ) * 0.5 / dr_f[ir-1];
            pBarInv_f = tmp_pBarInv_f + (pBarInv_f - tmp_pBarInv_f) * 0.5*dr[ir-1]/dr_f[ir-1];
        }
        // This is for uniform step sizes:
        // pBarInv_f += tmp_pBarInv_f;
        // pBarInv_f *= 0.5;
    }

    return pBarInv_f;
}




/**
 * Function for evalutaing a vector, defined in the cell grid, on the
 * flux grid. Interpolating in the inner, extrapolating on the outer
 * end, and setting zero flux on the inner most point.
 */ 
template<typename T>
real_t DREAM::SvenssonTransport<T>::EvalOnFluxGrid(len_t ir, const real_t *vec){
    // Need interpolation from cell grid to flux grid:
    // vec_f[0]  = vec[0]
    // vec_f[ir] = interpolate
    // vec_f[nr] = extrapolate


    // Inverse of p-bar on the Flux grid, with additional helper variable.
    real_t vec_f; 

    // Grid step size in the radial grid for the derivative.
    const real_t *dr_f = this->grid->GetRadialGrid()->GetDr_f();
    const real_t *dr = this->grid->GetRadialGrid()->GetDr(); 

    
    // Interpolating (extrapolating) the inverse of p bar onto the
    // flux grid.
    if (ir == 0) {
        // Zero flux at r = 0. Therefore choose the value at "ir=1/2".
        vec_f = vec[0];
    }
    else if (ir == this->nr_f - 1) {
        // Linearly extrapolating the value at the end point from the
        // two previous points.
        vec_f = vec[ir-1] + (vec[ir-1] - vec[ir-2]) * 0.5*dr[ir-1]/dr_f[ir-2];        
    }
    else {
        // In the middle, we simply linearly interpolate
        vec_f = vec[ir-1] + (vec[ir] - vec[ir-1]) * 0.5*dr[ir]/dr_f[ir-1];
    }

    return vec_f;
}
