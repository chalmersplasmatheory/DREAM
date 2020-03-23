#ifndef _DREAM_FVM_EQUATION_HPP
#define _DREAM_FVM_EQUATION_HPP

#include <vector>
#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAM::FVM {
    class Equation {
    private:
        std::vector<BC::BoundaryCondition*> boundaryConditions;
        std::vector<EquationTerm> terms;
        RadialGrid *grid;

    public:
        ~Equation();

        void RebuildTerms(const real_t);
        void SetMatrixElements(Matrix*);
    };
}

#endif/*_DREAM_FVM_EQUATION_HPP*/
