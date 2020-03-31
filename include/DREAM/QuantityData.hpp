#ifndef _DREAM_QUANTITY_DATA_HPP
#define _DREAM_QUANTITY_DATA_HPP

#include <petsc.h>
#include <vector>
#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"

namespace DREAM {
    class QuantityData {
    private:
        FVM::Grid *grid;
        std::vector<real_t> times;
        std::vector<real_t*> store;

        len_t nElements=0;
        real_t *data=nullptr;

        // Vector used for addressing PETSc vectors
        PetscInt *idxVec = nullptr;

        void AllocateData();

    public:
        QuantityData(FVM::Grid*);
        ~QuantityData();

        real_t *Get() { return this->data; }
        len_t Size() { return this->nElements; }

        void SaveStep(const real_t);
        void Store(const len_t, Vec&);
    };
}

#endif/*_DREAM_QUANTITY_DATA_HPP*/
