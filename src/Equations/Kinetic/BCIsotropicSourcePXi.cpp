/**
 * Implementation of the boundary condition used on the isotropic hot-tail
 * grid to facilitate transfer of particles from/to the cold electron
 * population.
 *
 * NOTE: This boundary condition only works for p/xi grids.
 */

#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/Kinetic/BCIsotropicSourcePXi.hpp"


using namespace DREAM;



/**
 * Constructor.
 */
BCIsotropicSourcePXi::BCIsotropicSourcePXi(FVM::Grid *g, CollisionQuantityHandler *cqh, len_t id_f)
    : FVM::BC::PInternalBoundaryCondition(g), id_f(id_f) {
    
    SetName("BCIsotropicSourcePXi");
    this->slowingDownFreq = cqh->GetNuS();
}


/**
 * Rebuild the flux.
 */
bool BCIsotropicSourcePXi::Rebuild(const real_t, FVM::UnknownQuantityHandler *uqh) {
    const real_t *f = uqh->GetUnknownData(this->id_f);

    len_t offset=0;
    for (len_t ir = 0; ir < grid->GetNr(); ir++) {
        const len_t np = grid->GetMomentumGrid(ir)->GetNp1();
        const len_t nxi = grid->GetMomentumGrid(ir)->GetNp2();
        const real_t p3nuS = this->slowingDownFreq->GetP3NuSAtZero(ir);
        const real_t *Vp_p2 = this->grid->GetVpOverP2AtZero(ir);

        for (len_t j = 0; j < nxi; j++)
            this->VpS[ir][j] = p3nuS*Vp_p2[j] * f[offset + np*j];
        offset += np * nxi;
    }

    return true;
}

