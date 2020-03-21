#ifndef _TQS_FVM_EQUATION_SYSTEM_HPP
#define _TQS_FVM_EQUATION_SYSTEM_HPP

namespace TQS::FVM { class EquationSystem; }

#include <petscmat.h>
#include <vector>
#include "FVM/config.h"
#include "FVM/FVMException.hpp"
#include "FVM/Matrix.hpp"

namespace TQS::FVM {
    class EquationSystem : public Matrix {
        private:
            struct _subeq {
                PetscInt offset;
                PetscInt n;     // Number of elements in unknown vector
                PetscInt nnz;   // Number of non-zero elements in corresponding part of matrix
            };
            std::vector<struct _subeq> subeqs;
            PetscInt next_subindex=0;

        public:
            EquationSystem();
            ~EquationSystem();

            // Block API
            void ConstructSystem();
            len_t CreateSubEquation(const PetscInt, const PetscInt);
            void SelectSubEquation(const PetscInt, const PetscInt);
            void RestoreSubEquation(Matrix*, const PetscInt, const PetscInt);

            void ZeroEquation(const PetscInt);
    };

    class EquationSystemException : public FVMException {
    public:
        template<typename ... Args>
        EquationSystemException(const std::string &msg, Args&& ... args)
            : MatrixException(msg, std::forward<Args>(args) ...) {
            AddModule("EquationSystem");
        }
    };
}

#endif/*_TQS_FVM_EQUATION_SYSTEM_HPP*/

