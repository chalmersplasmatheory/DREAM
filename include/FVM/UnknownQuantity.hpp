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

        // This variable can be used to represent multiple unknowns (i.e.
        // of size equal to the grid size) by a single unknown.
        len_t nMultiples=1;

    public:
        UnknownQuantity(const std::string& name, Grid *grid, const len_t nMultiples=1) {
            this->name = name;
            this->grid = grid;
            this->data = new QuantityData(grid, nMultiples);
            this->nMultiples = nMultiples;
        }
        ~UnknownQuantity() {
            delete data;
        }

        real_t *GetData() { return this->data->Get(); }
        real_t *GetDataPrevious() { return this->data->GetPrevious(); }
        real_t *GetInitialData() { return this->data->GetInitialData(); }
        const Grid *GetGrid() const { return this->grid; }
        const std::string& GetName() const { return this->name; }

        bool HasChanged() const { return data->HasChanged(); }
        bool HasInitialValue() const { return data->HasInitialValue(); }

        len_t NumberOfElements() const { return grid->GetNCells() * this->nMultiples; }

        void SaveStep(const real_t t, bool trueSave) { data->SaveStep(t, trueSave); }
        void Store(Vec& v, const len_t offs, bool mayBeConstant=false) { data->Store(v, offs, mayBeConstant); }
        void Store(const real_t *v, const len_t offs=0, bool mayBeConstant=false) { data->Store(v, offs, mayBeConstant); }

        void SaveSFile(SFile *sf, const std::string& path="", bool saveMeta=false)
        { this->data->SaveSFile(sf, this->name, path, saveMeta); }

        void SetInitialValue(const real_t*, const real_t t0=0);
    };
}

#endif/*_DREAM_FVM_UNKNOWN_QUANTITY_HPP*/
