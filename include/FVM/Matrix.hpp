#ifndef _DREAM_FVM_MATRIX_HPP
#define _DREAM_FVM_MATRIX_HPP

#include <petscis.h>
#include <petscmat.h>
#include <vector>
#include "FVM/config.h"
#include "FVM/FVMException.hpp"

namespace DREAM::FVM {
    class Matrix {
        protected:
            Mat petsc_mat;
            PetscInt m, n;
            PetscInt nz, *nnz;

            PetscInt rowOffset=0, colOffset=0;

            bool allocated=false;

            void Construct(
                const PetscInt, const PetscInt,
                const PetscInt, const PetscInt* nnzl=nullptr
            );

        public:
            enum view_format {
                ASCII_MATLAB      = PETSC_VIEWER_ASCII_MATLAB,
                ASCII_DENSE       = PETSC_VIEWER_ASCII_DENSE,
                ASCII_COMMON      = PETSC_VIEWER_ASCII_COMMON,
                ASCII_INFO        = PETSC_VIEWER_ASCII_INFO,
                ASCII_INFO_DETAIL = PETSC_VIEWER_ASCII_INFO_DETAIL,
                BINARY_MATLAB     = PETSC_VIEWER_BINARY_MATLAB,
                NON_ZERO_STRUCT   = 1337
            };

            Matrix();
            Matrix(const PetscInt, const PetscInt, Mat);
            Matrix(
                const PetscInt, const PetscInt,
                const PetscInt, const PetscInt *nnzl=nullptr
            );
            virtual ~Matrix();

            void Assemble();
            void PartialAssemble();
            bool ContainsNaNOrInf(len_t *I=nullptr, len_t *J=nullptr);
            void DiagonalScale(Vec, Vec);
            virtual void Destroy();
            void GetOwnershipRange(PetscInt*, PetscInt*);
            virtual void IMinusDtA(const PetscScalar);
            real_t *Multiply(const len_t, const real_t*);

            real_t GetElement(const PetscInt, const PetscInt);
            void GetRow(const PetscInt, PetscScalar*);
            void GetRow(const PetscInt, const PetscInt*, PetscScalar*);
            void GetColumn(const PetscInt, PetscScalar*);
            void GetColumn(const PetscInt, const PetscInt*, PetscScalar*);
            real_t *GetRowMaxAbs();
            void GetRowMaxAbs(real_t*);
            len_t GetNRows() const { return this->m; }
            len_t GetNCols() const { return this->n; }
            len_t GetNNZ();

            PetscInt GetRowOffset() const { return this->rowOffset; }
            PetscInt GetColOffset() const { return this->colOffset; }

            void SetElement(
                const PetscInt, const PetscInt,
                const PetscScalar, InsertMode im=ADD_VALUES
            );
			void SetRow(
				PetscInt, const PetscInt,
				PetscInt*, const PetscScalar*,
				InsertMode im=ADD_VALUES
			);

            void ResetOffset();
            void SetOffset(const PetscInt, const PetscInt);
            void View(enum view_format vf=ASCII_MATLAB, const std::string& filename="petsc_matrix");
            void Zero(bool nzKeep = true);
			void ZeroRows(const PetscInt, const PetscInt[]);
			void ZeroRowsColumns(const PetscInt, const PetscInt[]);
			void SetDiagonalConstant(const PetscInt, const PetscInt[], const PetscReal);

            void PrintInfo();

            Mat &mat() { return this->petsc_mat; }
    };

    class MatrixException : public FVMException {
    public:
        template<typename ... Args>
        MatrixException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("Matrix");
        }
    };
}

#endif/*_DREAM_FVM_MATRIX_HPP*/

