#ifndef _ADAS_DATA_H
#define _ADAS_DATA_H
/**
 * This file contains declarations of ADAS data and
 * data structures. This file is independent of DREAM.
 */

#include "FVM/config.h"

struct adas_rate {
    const char *name;
    len_t Z, A;

    // ADAS data
    len_t acd_nn, acd_nT;
    const real_t *acd_n, *acd_T, *acd;

    len_t ccd_nn, ccd_nT;
    const real_t *ccd_n, *ccd_T, *ccd;

    len_t scd_nn, scd_nT;
    const real_t *scd_n, *scd_T, *scd;

    len_t plt_nn, plt_nT;
    const real_t *plt_n, *plt_T, *plt;

    len_t prb_nn, prb_nT;
    const real_t *prb_n, *prb_T, *prb;
};

extern const len_t adas_rate_n;
extern struct adas_rate adas_rate_table[];


#endif/*_ADAS_DATA_H*/
