#ifndef _DREAM_OTHER_QUANTITY_HPP
#define _DREAM_OTHER_QUANTITY_HPP

#include <functional>
#include <string>
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/fluxGridType.enum.hpp"
#include "FVM/QuantityData.hpp"


namespace DREAM {
    class OtherQuantity {
    private:
        FVM::QuantityData *data;

        std::string name, description;
        FVM::Grid *grid;
        len_t nMultiples;
        enum FVM::fluxGridType fgt;

        bool active = false;

        std::function<void(FVM::QuantityData*)> storeFunc;

    public:
        OtherQuantity(
            const std::string& name, const std::string& desc,
            FVM::Grid *grid, const len_t nMultiples, enum FVM::fluxGridType fgt,
            std::function<void(FVM::QuantityData*)> storeFunc
        ) : name(name), description(desc), grid(grid), nMultiples(nMultiples), fgt(fgt) {

            this->storeFunc = storeFunc;
        }
        ~OtherQuantity() {
            delete data;
        }

        void Activate() { 
            this->data = new FVM::QuantityData(grid, nMultiples, fgt);  
            this->active = true;    
        }
        bool IsActive() { return this->active; }
        const std::string& GetName() { return this->name; }

        FVM::Grid *GetGrid() { return this->grid; }

        void SaveSFile(SFile *sf, const std::string& path="") {
            this->data->SaveSFile(sf, this->name, path, this->description);
        }
        void Store(const real_t t) {
            this->storeFunc(this->data);
            this->data->SaveStep(t, true);
        }
    };
}

#endif/*_DREAM_OTHER_QUANTITY_HPP*/
