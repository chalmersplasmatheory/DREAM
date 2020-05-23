/**
 * Implementation of the p/xi grid generator.
 */

#include <cmath>
#include "FVM/config.h"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"


using namespace DREAM::FVM::PXiGrid;


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
    const real_t t, const len_t ri, DREAM::FVM::MomentumGrid *mg,
    const DREAM::FVM::RadialGrid *rg
) {
    bool built = this->pGenerator->Rebuild(t, ri, mg, rg);
    built |= this->xiGenerator->Rebuild(t, ri, mg, rg);

    len_t np1 = mg->GetNp1();
    len_t np2 = mg->GetNp2();
    
    real_t
        *p    = new real_t[np1*np2],
        *p_f1 = new real_t[(np1+1)*np2],
        *p_f2 = new real_t[np1*(np2+1)],
        *gamma    = new real_t[np1*np2],
        *gamma_f1 = new real_t[(np1+1)*np2],
        *gamma_f2 = new real_t[np1*(np2+1)],
        *xi0  = new real_t[np1*np2],
        *xi01 = new real_t[(np1+1)*np2],
        *xi02 = new real_t[np1*(np2+1)];
    
    const real_t *p1   = mg->GetP1();
    const real_t *p1_f = mg->GetP1_f();
    const real_t *p2   = mg->GetP2();
    const real_t *p2_f = mg->GetP2_f();
    
    len_t pind;
    for (len_t j = 0; j < np2; j++) {
        for (len_t i = 0; i < np1; i++) {
            pind = np1*j+i;
            p[pind]   = p1[i];
            gamma[pind] = sqrt(1+p[pind]*p[pind]);
            xi0[pind] = p2[j];
        }
    }
    for (len_t j = 0; j < np2; j++) {
        for (len_t i = 0; i < np1+1; i++) {
            pind = (np1+1)*j+i;
            p_f1[pind] = p1_f[i];
            gamma_f1[pind] = sqrt(1+p_f1[pind]*p_f1[pind]);
            xi01[pind] = p2[j];
        }
    }
    for (len_t j = 0; j < np2+1; j++) {
        for (len_t i = 0; i < np1; i++) {
            pind = np1*j+i;
            p_f2[pind] = p1[i];
            gamma_f2[pind] = sqrt(1+p_f2[pind]*p_f2[pind]);
            xi02[pind] = p2_f[j];
        }
    }
    

    mg->InitializePAndXi0(p, p_f1, p_f2, gamma, gamma_f1, gamma_f2, xi0, xi01, xi02);

    return built;

}

