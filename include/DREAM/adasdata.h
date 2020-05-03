#ifndef _ADAS_DATA_H
#define _ADAS_DATA_H
/**
 * This file contains declarations of ADAS data and
 * data structures. This file is independent of DREAM.
 */

#include "FVM/config.h"

struct adas_rate {
    const char *name;
    len_t Z;

    // ADAS data
    len_t acd_nn, acd_nT;
    const real_t *acd, *acd_n, *acd_T;

    len_t scd_nn, scd_nT;
    const real_t *scd, *scd_n, *scd_T;

    len_t plt_nn, plt_nT;
    const real_t *plt, *plt_n, *plt_T;

    len_t prb_nn, prb_nT;
    const real_t *prb, *prb_n, *prb_T;
};

extern const len_t adas_rate_n;
extern struct adas_rate adas_rate_table[];


#endif/*_ADAS_DATA_H*/
