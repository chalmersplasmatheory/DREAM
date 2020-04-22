#ifndef _DREAM_QUANTITY_DATA_HPP
#define _DREAM_QUANTITY_DATA_HPP

#include <petsc.h>
#include <vector>
#include <softlib/SFile.h>
#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"

namespace DREAM::FVM {
    class QuantityData {
    private:
        Grid *grid;
        std::vector<real_t> times;
        std::vector<real_t*> store;

        len_t nElements=0;
        real_t *data=nullptr;

        // Vector used for addressing PETSc vectors
        PetscInt *idxVec = nullptr;

        bool hasChanged = true;

        void AllocateData();

    public:
        QuantityData(FVM::Grid*);
        ~QuantityData();

        real_t *Get() { return this->data; }
        real_t *GetPrevious() { return this->store.back(); }
        len_t Size() { return this->nElements; }

        /**
         * Returns 'true' if the data stored by this quantity
         * was changed in the last call to 'Store()'.
         */
        bool HasChanged() const { return this->hasChanged; }
        bool HasInitialValue() const { return (this->store.size()>=1); }

        void SaveStep(const real_t);
        void Store(Vec&, const len_t, bool mayBeConstant=false);
        void Store(const real_t*, const len_t offset=0, bool mayBeConstant=false);

        void SaveSFile(SFile*, const std::string& name, const std::string& path="", bool saveMeta=false);

        void SetInitialValue(const real_t*, const real_t t0=0);
    };
}

#endif/*_DREAM_QUANTITY_DATA_HPP*/
