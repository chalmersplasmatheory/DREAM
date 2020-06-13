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

extern const len_t nist_binding_n;
extern struct nist_data nist_binding_table[];

extern const len_t nist_ionization_n;
extern struct nist_data nist_ionization_table[];

#endif/*_NIST_DATA_H*/
