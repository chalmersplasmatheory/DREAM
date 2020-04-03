#ifndef _DREAM_FVM_EQUATION_HPP
#define _DREAM_FVM_EQUATION_HPP

#include <vector>
#include "FVM/Equation/AdvectionDiffusionTerm.hpp"
#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Grid/Grid.hpp"

namespace DREAM::FVM {
    class EquationException : public FVMException {
    public:
        template<typename ... Args>
        EquationException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("Equation");
        }
    };

    class Equation {
    private:
        std::vector<BC::BoundaryCondition*> boundaryConditions;
        std::vector<EquationTerm*> terms;
        AdvectionDiffusionTerm *adterm = nullptr;
        TransientTerm *tterm = nullptr;
        Grid *grid;

        enum AdvectionDiffusionTerm::advdiff_interpolation advdiff_interpolationMethod;

    public:
        Equation(
            Grid*, enum AdvectionDiffusionTerm::advdiff_interpolation ip=AdvectionDiffusionTerm::AD_INTERP_CENTRED
        );

        ~Equation();

        void AddTerm(AdvectionTerm *a) {
            if (adterm == nullptr)
                adterm = new AdvectionDiffusionTerm(this->grid, this->advdiff_interpolationMethod);

            adterm->Add(a);
        }
        void AddTerm(DiffusionTerm *d) {
            if (adterm == nullptr)
                adterm = new AdvectionDiffusionTerm(this->grid, this->advdiff_interpolationMethod);

            adterm->Add(d);
        }
        void AddTerm(TransientTerm *t) {
            if (tterm != nullptr)
                throw EquationException("The equation already has a transient term.");

            tterm = t;
        }
        void AddTerm(EquationTerm *t)  { terms.push_back(t); }

        void AddBoundaryCondition(BC::BoundaryCondition *bc) {
            boundaryConditions.push_back(bc);
        }

        void RebuildTerms(const real_t, const real_t);
        void SetMatrixElements(Matrix*, real_t*);
    };
}

#endif/*_DREAM_FVM_EQUATION_HPP*/
