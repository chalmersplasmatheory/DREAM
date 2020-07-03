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
                PetscInt id;    // Externally set ID which can be used to identify the block
                PetscInt offset;// Row/col offset of block in matrix
                PetscInt n;     // Number of elements in unknown vector
                PetscInt nnz;   // Number of non-zero elements in corresponding part of matrix
            };
            std::vector<struct _subeq> subeqs;
            PetscInt next_subindex=0;

            PetscInt blockn=0;

        public:
            BlockMatrix();
            virtual ~BlockMatrix();

            // Block API
            void ConstructSystem();
            len_t CreateSubEquation(const PetscInt, const PetscInt, const PetscInt id=-1);
            PetscInt GetOffset(const PetscInt);
            PetscInt GetOffsetById(const PetscInt);
            void SelectSubEquation(const PetscInt, const PetscInt);
            void RestoreSubEquation(Matrix*, const PetscInt, const PetscInt);

            len_t GetNNZInBlock(const len_t i) const { return subeqs.at(i).nnz; }

            virtual void IMinusDtA(const PetscScalar) override;

            void ZeroEquation(const PetscInt);
    };

    class BlockMatrixException : public MatrixException {
    public:
        template<typename ... Args>
        BlockMatrixException(const std::string &msg, Args&& ... args)
            : MatrixException(msg, std::forward<Args>(args) ...) {
            AddModule("BlockMatrix");
        }
    };
}

#endif/*_DREAM_FVM_BLOCK_MATRIX_HPP*/

