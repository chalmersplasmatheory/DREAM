#ifndef _DREAM_QUANTITY_DATA_HPP
#define _DREAM_QUANTITY_DATA_HPP

#include <petsc.h>
#include <vector>
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

        void AllocateData();

    public:
        QuantityData(FVM::Grid*);
        ~QuantityData();

        real_t *Get() { return this->data; }
        len_t Size() { return this->nElements; }

        void SaveStep(const real_t);
        void Store(Vec&, const len_t);
        void Store(const real_t*, const len_t offset=0);
    };
}

#endif/*_DREAM_QUANTITY_DATA_HPP*/
