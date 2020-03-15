#ifndef _TQS_FVM_EQUATION_HPP
#define _TQS_FVM_EQUATION_HPP

#include <vector>
#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace TQS::FVM {
    class Equation {
    private:
        std::vector<BoundaryCondition> boundaryConditions;
        std::vector<EquationTerm> terms;
        RadialGrid *grid;

    public:
        ~Equation();

        void RebuildTerms(const real_t);
        void SetMatrixElements(Matrix*);
    };
}

#endif/*_TQS_FVM_EQUATION_HPP*/
