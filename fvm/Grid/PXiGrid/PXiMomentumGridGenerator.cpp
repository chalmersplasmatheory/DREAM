/**
 * Implementation of the p/xi grid generator.
 */

#include <cmath>
#include "config.h"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"


using namespace TQS::FVM::PXiGrid;


/**
 * Destructor.
 */
MomentumGridGenerator::~MomentumGridGenerator() {
    delete this->pGenerator;
    delete this->xiGenerator;
}

/**
 * Returns true if this momentum grid needs to
 * be re-built.
 *
 * t:            For which the grid should be re-built.
 * rGridRebuilt: True if the radial grid associated with
 *               with this generator was re-built for this
 *               time step.
 */
bool MomentumGridGenerator::NeedsRebuild(
    const real_t t, const bool rGridRebuilt
) {
    return (
        this->pGenerator->NeedsRebuild(t, rGridRebuilt) ||
        this->xiGenerator->NeedsRebuild(t, rGridRebuilt)
    );
}

/**
 * Re-builds the given momentum grid using this
 * momentum grid generator.
 */
bool MomentumGridGenerator::Rebuild(
    const real_t t, const len_t ri, TQS::FVM::MomentumGrid *mg,
    const TQS::FVM::RadialGrid *rg
) {
    bool built = this->pGenerator->Rebuild(t, ri, mg, rg);
    built |= this->xiGenerator->Rebuild(t, ri, mg, rg);

    if (!built)
        return false;

    const len_t np  = this->pGenerator->GetNp();
    const len_t nxi = this->xiGenerator->GetNxi();
    const len_t N  = np*nxi;

    const real_t
        *p    = mg->GetP1(),
        *p_f  = mg->GetP1_f(),
        *xi   = mg->GetP2(),
        *xi_f = mg->GetP2_f(),
        *dp   = mg->GetDp1(),
        *dxi  = mg->GetDp2();

    real_t
        *volumes = new real_t[N];

    // Construct grid volumes
    for (len_t j = 0; j < nxi; j++) {
        for (len_t i = 0; i < np; i++)
            volumes[j*np + i] = 2*M_PI*p[i]*p[i] * dp[i] * dxi[j];
    }

    // Construct scale factors
    real_t
        *hp, *hxi, *hphi,
        *hp_fp, *hxi_fp, *hphi_fp,
        *hp_fx, *hxi_fx, *hphi_fx;

    GenerateLameCoeffs(np, nxi, p, xi, &hp, &hxi, &hphi);
    GenerateLameCoeffs(np+1, nxi, p_f, xi, &hp_fp, &hxi_fp, &hphi_fp);
    GenerateLameCoeffs(np, nxi+1, p, xi_f, &hp_fx, &hxi_fx, &hphi_fx);

    mg->InitializeMetric(
        volumes,
        hp, hxi, hphi,
        hp_fp, hxi_fp, hphi_fp,
        hp_fx, hxi_fx, hphi_fx
    );

    return built;
}

/**
 * Generate Lamé coefficients on the given grid.
 *
 * np:   Number of points in momentum grid.
 * nxi:  Number of points in pitch grid.
 * p:    Momentum grid.
 * xi:   Pitch grid.
 * hp:   Lamé coefficient for momentum coordinate.
 * hxi:  Lamé coefficient for pitch coordinate.
 * hphi: Lamé coefficient for gyro phase coordinate.
 */
void MomentumGridGenerator::GenerateLameCoeffs(
    const len_t np, const len_t nxi,
    const real_t *p, const real_t *xi,
    real_t **hp, real_t **hxi, real_t **hphi
) {
    *hp   = new real_t[np*nxi];
    *hxi  = new real_t[np*nxi];
    *hphi = new real_t[np*nxi];

    for (len_t j = 0; j < nxi; j++) {
        for (len_t i = 0; i < np; i++) {
            real_t sint = sqrt(1.0 - xi[j]*xi[j]);

            (*hp)  [j*np + i] = 1;
            (*hxi) [j*np + i] =-p[i] / sint;
            (*hphi)[j*np + i] = p[i] * sint;
        }
    }
}

