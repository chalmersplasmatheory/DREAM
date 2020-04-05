#ifndef _DREAM_FVM_UNKNOWN_QUANTITY_HANDLER_HPP
#define _DREAM_FVM_UNKNOWN_QUANTITY_HANDLER_HPP

namespace DREAM::FVM { class UnknownQuantityHandler; }

#include <string>
#include <vector>
#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantity.hpp"

namespace DREAM::FVM {
    class UnknownQuantityHandler {
    private:
        std::vector<UnknownQuantity*> unknowns;

    public:
        UnknownQuantityHandler();
        ~UnknownQuantityHandler();

        UnknownQuantity *GetUnknown(const len_t i) { return unknowns.at(i); }
        len_t GetUnknownID(const std::string&);
        len_t GetNUnknowns() const { return this->unknowns.size(); }
        len_t Size() const { return GetNUnknowns(); }

        real_t *GetUnknownData(const len_t);
        real_t *GetUnknownDataPrevious(const len_t);

        len_t InsertUnknown(const std::string&, Grid*);

        UnknownQuantity *operator[](const len_t i) { return GetUnknown(i); }
    };
}

#endif/*_DREAM_FVM_UNKNOWN_QUANTITY_HANDLER_HPP*/
