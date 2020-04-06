#ifndef _DREAM_FVM_UNKNOWN_QUANTITY_HPP
#define _DREAM_FVM_UNKNOWN_QUANTITY_HPP

#include <map>
#include <string>
#include "FVM/Grid/Grid.hpp"
#include "FVM/QuantityData.hpp"

namespace DREAM::FVM {
    class UnknownQuantity {
    private:
        // Name of quantity
        std::string name;
        // Pointer to grid on which the quantity is defined
        Grid *grid;
        // (Solution) data handler
        QuantityData *data;

    public:
        UnknownQuantity(const std::string& name, Grid *grid) {
            this->name = name;
            this->grid = grid;
            this->data = new QuantityData(grid);
        }
        ~UnknownQuantity() {
            delete data;
        }

        real_t *GetData() { return this->data->Get(); }
        const Grid *GetGrid() const { return this->grid; }
        const std::string& GetName() const { return this->name; }

        len_t NumberOfElements() const { return grid->GetNCells(); }

        void SaveStep(const real_t t) { data->SaveStep(t); }
        void Store(Vec& v, const len_t offs) { data->Store(v, offs); }
        void Store(const real_t *v, const len_t offs=0) { data->Store(v, offs); }
    };
}

#endif/*_DREAM_FVM_UNKNOWN_QUANTITY_HPP*/
