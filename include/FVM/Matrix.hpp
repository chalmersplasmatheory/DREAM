#ifndef _DREAM_FVM_MATRIX_HPP
#define _DREAM_FVM_MATRIX_HPP

#include <petscis.h>
#include <petscmat.h>
#include <vector>
#include "FVM/config.h"
#include "FVM/FVMException.hpp"

namespace DREAM::FVM {
    /**
     * Matrix wrapper around a PETSc Mat.
     *
     * The Matrix object may represent a view into a sub-block of a larger
     * matrix. In that case, rowOffset/colOffset are added to all indices
     * passed to setters/getters.
     *
     * Performance:
     *  - Repeated calls to MatSetValue() can be expensive. For assembly code
     *    that inserts many entries per row, the buffered insertion API can
     *    reduce PETSc call overhead by batching row insertions via MatSetValues().
     */
    class Matrix {
        protected:
            Mat petsc_mat;
            PetscInt m, n;
            PetscInt nz, *nnz;

            PetscInt rowOffset=0, colOffset=0;

            std::vector<PetscInt> tmpRows;
            std::vector<PetscScalar> tmpVals;

            bool allocated=false;

            bool buffering = false;
            PetscInt bufRow = -1;
            PetscInt bufLen = 0;
            InsertMode bufMode = ADD_VALUES;
            std::vector<PetscInt> bufCols;
            std::vector<PetscScalar> bufVals;
            void FlushBufferedRow();

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

            /**
             * Set one matrix element.
             *
             * If buffered insertion is enabled (BeginBufferedSet), values are queued
             * and flushed per row using MatSetValues(). Otherwise, this calls MatSetValue().
             */
            void SetElement(
                const PetscInt, const PetscInt,
                const PetscScalar, InsertMode im=ADD_VALUES
            );

            /**
             * Insert values in a single column for a contiguous block of rows.
             *
             * This is a fast path for column-oriented assembly where the caller has a
             * contiguous block of values corresponding to rows:
             *
             *   row = firstLocalRow + r,   r = 0..m-1
             *   col = localCol
             *   value(row, col) = vals[r]
             *
             * Applies rowOffset/colOffset internally. Exact zeros are skipped to avoid
             * creating structural nonzeros for inactive rows (matching SetElement()).
             */
            void SetColumn(
                PetscInt firstLocalRow, PetscInt localCol,
                PetscInt m, const PetscScalar *vals,
                InsertMode im=ADD_VALUES
            );

            /**
             * Enable buffered insertion mode.
             *
             * While buffering is enabled, SetElement() accumulates (col,value) pairs
             * for the current row and flushes them with MatSetValues() when the row
             * changes or EndBufferedSet() is called.
             *
             * mode:    Insert mode for the whole buffered region (typically ADD_VALUES).
             * reserve: Optional hint for expected nnz per row.
             */
            void BeginBufferedSet(InsertMode mode = ADD_VALUES, PetscInt reserve = 0);

            /**
             * Disable buffered insertion mode and flush any pending row.
             */
            void EndBufferedSet();

            /**
             * Insert values in a single row for multiple columns (row-oriented bulk insert).
             *
             * Notes:
             *  - Applies rowOffset/colOffset.
             *  - Does not filter zeros (callers may intentionally insert zeros to
             *    establish a stable sparsity pattern).
             */
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

