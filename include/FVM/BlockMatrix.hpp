#ifndef _DREAM_FVM_BLOCK_MATRIX_HPP
#define _DREAM_FVM_BLOCK_MATRIX_HPP

namespace DREAM::FVM { class BlockMatrix; }

#include <petscmat.h>
#include <vector>
#include "FVM/config.h"
#include "FVM/FVMException.hpp"
#include "FVM/Matrix.hpp"

namespace DREAM::FVM {
    class BlockMatrix : public Matrix {
        private:
            struct _subeq {
                PetscInt offset;
                PetscInt n;     // Number of elements in unknown vector
                PetscInt nnz;   // Number of non-zero elements in corresponding part of matrix
            };
            std::vector<struct _subeq> subeqs;
            PetscInt next_subindex=0;

            PetscInt blockn=0;

        public:
            BlockMatrix();
            ~BlockMatrix();

            // Block API
            void ConstructSystem();
            len_t CreateSubEquation(const PetscInt, const PetscInt);
            void SelectSubEquation(const PetscInt, const PetscInt);
            void RestoreSubEquation(Matrix*, const PetscInt, const PetscInt);

            virtual void IMinusDtA(const PetscScalar) override;

            void ZeroEquation(const PetscInt);
    };

    class BlockMatrixException : public FVMException {
    public:
        template<typename ... Args>
        BlockMatrixException(const std::string &msg, Args&& ... args)
            : MatrixException(msg, std::forward<Args>(args) ...) {
            AddModule("BlockMatrix");
        }
    };
}

#endif/*_DREAM_FVM_BLOCK_MATRIX_HPP*/

