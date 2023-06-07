/**
 * Implementation of the hyperresistive diffusion of poloidal flux.
 * Based on the section on hyperresistivity in DREAM/doc/notes/theory
 */

#include "DREAM/Equations/Fluid/HyperresistiveDiffusionTerm.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
HyperresistiveDiffusionTerm::HyperresistiveDiffusionTerm(
    FVM::Grid *g, FVM::Interpolator1D *Lambda
) : FVM::DiffusionTerm(g), Lambda(Lambda) {
    
    SetName("HyperresistiveDiffusionTerm");
}

/**
 * Build the coefficients of this diffusion term.
 */
void HyperresistiveDiffusionTerm::Rebuild(const real_t t, const real_t, FVM::UnknownQuantityHandler *){
    FVM::RadialGrid *rGrid = grid->GetRadialGrid(); 
    const real_t *Lmbd  = this->Lambda->Eval(t);

    // XXX: here we assume that all radii have the same momentum grids
    const len_t np1 = n1[0], np2 = n2[0];

    // (skip ir=0 since psi_t=0 there, and to avoid 1/psitPrime = 1/0)
    for (len_t ir = 1; ir < nr+1; ir++) {
        real_t Bmin = rGrid->GetBmin_f(ir);
        real_t BdotPhi = rGrid->GetBTorG_f(ir)*rGrid->GetFSA_1OverR2_f(ir);
        real_t VpVol = rGrid->GetVpVol_f(ir); 

        real_t psitPrime = VpVol*BdotPhi / (2*M_PI);

        // The entire d psi/dt equation is multiplied by 2*pi*psi_t'/VpVol,
        // so we get an extra factor of 2*pi here, and only one factor
        // of psi_t' (the factor 1/VpVol from the normalisation is included 
        // via the DiffusionTerm class). We also add a factor of 1/VpVol 
        // to the diffusion coefficient to cancel the factor VpVol inside 
        // the first radial derivative added by the DiffusionTerm. 
        //
        // Also, we divide by 'Bmin' since this operator is applied to
        // 'j_tot / (B/Bmin)'.
        real_t drr = 
            2*M_PI*rGrid->GetToroidalFlux_f(ir)*Lmbd[ir] / (VpVol*psitPrime*Bmin);

        for (len_t j = 0; j < np2; j++) 
            for (len_t i = 0; i < np1; i++) 
                Drr(ir, i, j) += drr;
    }
}

