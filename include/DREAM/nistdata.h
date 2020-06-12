#ifndef _NIST_DATA_H
#define _NIST_DATA_H
/**
 * This file contains declarations of NIST ADS data and
 * data structures. This file is independent of DREAM.
 */

#include "FVM/config.h"

struct nist_data {
    const char *name;
    len_t Z;
    const real_t *data;
};

#endif/*_NIST_DATA_H*/
