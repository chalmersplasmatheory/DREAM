/**
 * Implementation of the hyperresistive diffusion of poloidal flux.
 * Based on the section on hyperresistivity in DREAM/doc/notes/theory
 */

#include "DREAM/Equations/Fluid/HyperresistiveDiffusionTerm.hpp"


using namespace DREAM;

/**
 * Constructor.
 * TODO: the constructor should probably take a "TransportCoefficientHandler", which has a method GetLambda or similar that Rebuild can call.
 *       psi_t should probably be retrieved from radialgrid  
 */
HyperresistiveDiffusionTerm::HyperresistiveDiffusionTerm(FVM::Grid *g, real_t *Lambda, real_t *psi_t) : 
    FVM::DiffusionTerm(g), Lambda(Lambda), psi_t(psi_t) { }

/**
 * Build the coefficients of this diffusion term.
 */
void HyperresistiveDiffusionTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *){
    FVM::RadialGrid *rGrid = grid->GetRadialGrid(); 
    for (len_t ir = 0; ir < nr; ir++) {
        real_t BdotPhi = rGrid->GetBTorG(ir)*rGrid->GetFSA_1OverR2(ir);
        real_t Bmin = rGrid->GetBmin(ir);
        real_t VpVol = rGrid->GetVpVol(ir); 
        real_t drr = 4*M_PI*M_PI*rGrid->GetToroidalFlux(ir)*Lambda[ir]/(VpVol*BdotPhi*Bmin*Bmin);
        for (len_t j = 0; j < n2[ir]; j++) 
            for (len_t i = 0; i < n1[ir]; i++) 
                Drr(ir, i, j) += drr;
    }
}

