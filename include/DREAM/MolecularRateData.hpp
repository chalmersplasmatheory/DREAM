#ifndef _DREAM__MOLECULARRATEDATA_HPP
  #define _DREAM__MOLECULARRATEDATA_HPP

  #include "FVM/config.h"

  namespace DREAM {
      struct molecular_rate {
          const char *name;

          len_t nn;
          len_t nT;

          const real_t *logn;
          const real_t *logT;
          const real_t *coeff;
      };

      extern const len_t molecular_rate_n;
      extern struct molecular_rate molecular_rate_table[];
  }

  #endif