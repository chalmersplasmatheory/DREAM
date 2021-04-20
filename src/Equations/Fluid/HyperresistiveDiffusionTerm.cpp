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

    for (len_t ir = 0; ir < nr+1; ir++) {
        real_t Bmin = rGrid->GetBmin_f(ir);
        real_t BdotPhi = rGrid->GetBTorG_f(ir)*rGrid->GetFSA_1OverR2_f(ir);
        real_t VpVol = rGrid->GetVpVol_f(ir); 

        real_t psitPrime = VpVol*BdotPhi / (2*M_PI);

        // The entire d psi/dt equation is multiplied by 2*pi*psi_t',
        // so we get an extra factor of 2*pi here, and only one factor
        // of psi_t'.
        //
        // Also, we divide by 'Bmin' since this operator is applied to
        // 'j_tot / (B/Bmin)'.
        real_t drr = 
            2*M_PI*rGrid->GetToroidalFlux_f(ir)*Lmbd[ir] / (psitPrime*Bmin);

        for (len_t j = 0; j < n2[ir]; j++) 
            for (len_t i = 0; i < n1[ir]; i++) 
                Drr(ir, i, j) += drr;
    }
}

