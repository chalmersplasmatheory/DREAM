#ifndef _DREAM_MOLLER_PITCH_CONTRACTION_HPP
#define _DREAM_MOLLER_PITCH_CONTRACTION_HPP

#include <algorithm>
#include "FVM/config.h"

namespace DREAM {
    struct MollerPitchInterpolation {
        len_t m0;
        real_t w1;
    };

    /**
     * Accumulate one outgoing momentum row. Regroup the EXISTING two-point
     * angular interpolation before multiplying by its dense pitch matrices:
     *   U_ml = sum_k sigma_k b_km W_kl; C_j += sum_ml D_mjl U_ml.
     * This only changes summation order. No numerical sparsity threshold is used.
     *
     * W[k*nxiP+l], D[m][l*nxiK+j], scratch[m*nxiP+l]. Scratch is caller-owned,
     * size nSlices*nxiP, and is reset here. C is accumulated, not reset.
     * Interpolation indices/weights come from the validated angular kernel.
     */
    inline void AccumulateMollerPitch(
        len_t npP, len_t nxiP, len_t nxiK, len_t nSlices,
        const MollerPitchInterpolation *interp, const real_t *sigma,
        const real_t *W, const real_t *const *D, real_t *scratch, real_t *C
    ) {
        std::fill(scratch,scratch+nSlices*nxiP,0);
        len_t first=nSlices, last=0;
        for (len_t k=0; k<npP; ++k) {
            if (sigma[k]==0) continue;
            const auto& t=interp[k];
            first=std::min(first,t.m0);
            last=std::max(last,t.m0+1);
            const real_t a0=sigma[k]*(1-t.w1), a1=sigma[k]*t.w1;
            const real_t *w=W+k*nxiP;
            real_t *u0=scratch+t.m0*nxiP, *u1=u0+nxiP;
            for (len_t l=0; l<nxiP; ++l) {
                u0[l]+=a0*w[l];
                u1[l]+=a1*w[l];
            }
        }
        if (first==nSlices) return;
        for (len_t m=first; m<=last; ++m) {
            const real_t *u=scratch+m*nxiP;
            for (len_t l=0; l<nxiP; ++l) {
                if (u[l]==0) continue;
                const real_t *column=D[m]+l*nxiK;
                for (len_t j=0; j<nxiK; ++j)
                    C[j]+=u[l]*column[j];
            }
        }
    }
}
#endif
