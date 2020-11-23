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
        // Description of what the quantity represents
        std::string description;
        // Description of equation used to solve for this quantity
        std::string description_eqn;
        // Pointer to grid on which the quantity is defined
        Grid *grid;
        // (Solution) data handler
        QuantityData *data;

        // This variable can be used to represent multiple unknowns (i.e.
        // of size equal to the grid size) by a single unknown.
        len_t nMultiples=1;

    public:
        UnknownQuantity(
            const std::string& name, const std::string& description,
            Grid *grid, const len_t nMultiples=1
        ) {
            this->name = name;
            this->description = description;
            this->grid = grid;
            this->data = new QuantityData(grid, nMultiples);
            this->nMultiples = nMultiples;
        }
        ~UnknownQuantity() {
            delete data;
        }

        real_t GetPreviousTime() { return this->data->GetPreviousTime(); }
        real_t *GetData() { return this->data->Get(); }
        real_t *GetDataPrevious() { return this->data->GetPrevious(); }
        real_t *GetInitialData() { return this->data->GetInitialData(); }
        Grid *GetGrid() { return this->grid; }
        const std::string& GetDescription() const { return this->description; }
        const std::string& GetEquationDescription() const { return this->description_eqn; }
        const std::string& GetName() const { return this->name; }

        bool HasChanged() const { return data->HasChanged(); }
        bool HasInitialValue() const { return data->HasInitialValue(); }

        len_t NumberOfElements() const { return grid->GetNCells() * this->nMultiples; }
        len_t NumberOfMultiples() const { return this->nMultiples; }

        bool CanRollbackSaveStep() const { return data->CanRollbackSaveStep(); }
        void RollbackSaveStep() { data->RollbackSaveStep(); }
        void SaveStep(const real_t t, bool trueSave) { data->SaveStep(t, trueSave); }
        void SetEquationDescription(const std::string& d) { this->description_eqn = d; }
        void Store(Vec& v, const len_t offs, bool mayBeConstant=false) { data->Store(v, offs, mayBeConstant); }
        void Store(const real_t *v, const len_t offs=0, bool mayBeConstant=false) { data->Store(v, offs, mayBeConstant); }

        void SaveSFile(SFile *sf, const std::string& path="", bool saveMeta=false);
        void SaveSFileCurrent(SFile *sf, const std::string& path="", bool saveMeta=false);

        void SetInitialValue(const real_t*, const real_t t0=0);
    };
}

#endif/*_DREAM_FVM_UNKNOWN_QUANTITY_HPP*/
