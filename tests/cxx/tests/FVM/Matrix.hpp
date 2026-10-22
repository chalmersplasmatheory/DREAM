#ifndef _DREAMTESTS_FVM_MATRIX_HPP
#define _DREAMTESTS_FVM_MATRIX_HPP

#include <string>
#include <vector>
#include "FVM/Matrix.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::FVM {
    class Matrix : public UnitTest {
    public:
        struct entry {
            PetscInt row, col;
            PetscScalar val;
        };

        Matrix(const std::string& name) : UnitTest(name) {}

        bool CompareMatrices(DREAM::FVM::Matrix*, DREAM::FVM::Matrix*, const std::string&);
        bool HasStructuralEntry(DREAM::FVM::Matrix*, const PetscInt, const PetscInt);
        DREAM::FVM::Matrix *NewMatrix();
        std::vector<struct entry> GenerateEntries(bool);

        bool CheckBufferedAdd();
        bool CheckBufferedInsert();
        bool CheckBufferedFinalRowFlush();
        bool CheckSetColumn();

        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_FVM_MATRIX_HPP*/
