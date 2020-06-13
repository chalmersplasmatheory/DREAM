#ifndef _DREAM_NIST_HPP
#define _DREAM_NIST_HPP

#include <unordered_map>
#include "DREAM/nistdata.h"
#include "FVM/config.h"

namespace DREAM {
    class NIST {
    private:
        std::unordered_map<len_t, struct nist_data*> binding;
        std::unordered_map<len_t, struct nist_data*> ionization;

        real_t _GetEnergy(const len_t, const len_t, const std::unordered_map<len_t, struct nist_data*>&) const;
    public:
        NIST();

        real_t GetBindingEnergy(const len_t Z, const len_t Z0) const;
        real_t GetIonizationEnergy(const len_t Z, const len_t Z0) const;
    };
}

#endif/*_DREAM_NIST_HPP*/
