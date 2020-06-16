/**
 * Implementation of the diffusion term appearing in Ampere's law connecting plasma current and poloidal flux.
 */

#include "DREAM/Equations/PoloidalFlux/AmperesLawDiffusionTerm.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
AmperesLawDiffusionTerm::AmperesLawDiffusionTerm(FVM::Grid *g) : FVM::DiffusionTerm(g) { }

/**
 * Build the coefficients of this diffusion term.
 * XXX: assume same momentum grid at all radii
 */
void AmperesLawDiffusionTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *){
    for (len_t ir = 0; ir < nr+1; ir++) {
        real_t drr = grid->GetRadialGrid()->GetFSA_NablaR2OverR2_f(ir);
        for (len_t j = 0; j < n2[0]; j++) 
            for (len_t i = 0; i < n1[0]; i++) 
                Drr(ir, i, j) += drr;
    }
}

