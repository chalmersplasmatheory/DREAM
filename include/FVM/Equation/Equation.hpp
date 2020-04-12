#ifndef _DREAM_FVM_EQUATION_HPP
#define _DREAM_FVM_EQUATION_HPP

#include <vector>
#include "FVM/Equation/AdvectionDiffusionTerm.hpp"
#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
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
        PredeterminedParameter *predetermined = nullptr;
        AdvectionDiffusionTerm *adterm = nullptr;
        //TransientTerm *tterm = nullptr;
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

            CheckConsistency();
        }
        void AddTerm(DiffusionTerm *d) {
            if (adterm == nullptr)
                adterm = new AdvectionDiffusionTerm(this->grid, this->advdiff_interpolationMethod);

            adterm->Add(d);

            CheckConsistency();
        }
        void AddTerm(PredeterminedParameter *p) {
            if (predetermined != nullptr)
                throw EquationException("A predetermined parameter has already been applied to this quantity.");
            predetermined = p;

            CheckConsistency();
        }
        /*void AddTerm(TransientTerm *t) {
            if (tterm != nullptr)
                throw EquationException("The equation already has a transient term.");

            tterm = t;

            CheckConsistency();
        }*/
        void AddTerm(EquationTerm *t)  {
            terms.push_back(t);

            CheckConsistency();
        }

        void AddBoundaryCondition(BC::BoundaryCondition *bc) {
            boundaryConditions.push_back(bc);
        }

        // Verifies that the equation is consistent
        void CheckConsistency() {
            if (predetermined != nullptr) {
                if (adterm != nullptr/* || tterm != nullptr*/ || terms.size() > 0)
                    throw EquationException("A predetermined quantity cannot have other equation terms.");
            }
        }

        len_t GetNumberOfNonZerosPerRow() const;
        len_t GetNumberOfNonZerosPerRow_jac() const;
        PredeterminedParameter *GetPredetermined() { return this->predetermined; }
        /**
         * Returns 'true' if all terms of this equation are
         * 'PredeterminedParameter's.
         */
        bool IsPredetermined() const { return (predetermined != nullptr); }

        void RebuildTerms(const real_t, const real_t, UnknownQuantityHandler*);

        void SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix*);
        void SetMatrixElements(Matrix*, real_t*);
        void SetVectorElements(real_t*, const real_t*);
    };
}

#endif/*_DREAM_FVM_EQUATION_HPP*/
