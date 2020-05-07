#ifndef _DREAM_FVM_EQUATION_HPP
#define _DREAM_FVM_EQUATION_HPP

#include <vector>
#include "FVM/Equation/AdvectionDiffusionTerm.hpp"
#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/EvaluableEquationTerm.hpp"
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
        std::vector<EvaluableEquationTerm*> eval_terms;
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
        void AddTerm(EvaluableEquationTerm *t)  {
            eval_terms.push_back(t);

            CheckConsistency();
        }
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
                if (adterm != nullptr/* || tterm != nullptr*/ || terms.size() > 0 || boundaryConditions.size() > 0)
                    throw EquationException("A predetermined quantity cannot have other equation terms.");
            }
        }

        void Evaluate(real_t*, const real_t*);

        const real_t *const* GetAdvectionCoeffR() const { return this->adterm->GetAdvectionCoeffR(); }
        const real_t *GetAdvectionCoeffR(const len_t i) const { return this->adterm->GetAdvectionCoeffR(i); }
        const real_t *const* GetAdvectionCoeff1() const { return this->adterm->GetAdvectionCoeff1(); }
        const real_t *GetAdvectionCoeff1(const len_t i) const { return this->adterm->GetAdvectionCoeff1(i); }
        const real_t *const* GetAdvectionCoeff2() const { return this->adterm->GetAdvectionCoeff2(); }
        const real_t *GetAdvectionCoeff2(const len_t i) const { return this->adterm->GetAdvectionCoeff2(i); }

        const real_t *const* GetDiffusionCoeffRR() const { return this->adterm->GetDiffusionCoeffRR(); }
        const real_t *GetDiffusionCoeffRR(const len_t i) const { return this->adterm->GetDiffusionCoeffRR(i); }
        const real_t *const* GetDiffusionCoeff11() const { return this->adterm->GetDiffusionCoeff11(); }
        const real_t *GetDiffusionCoeff11(const len_t i) const { return this->adterm->GetDiffusionCoeff11(i); }
        const real_t *const* GetDiffusionCoeff12() const { return this->adterm->GetDiffusionCoeff12(); }
        const real_t *GetDiffusionCoeff12(const len_t i) const { return this->adterm->GetDiffusionCoeff12(i); }
        const real_t *const* GetDiffusionCoeff21() const { return this->adterm->GetDiffusionCoeff21(); }
        const real_t *GetDiffusionCoeff21(const len_t i) const { return this->adterm->GetDiffusionCoeff21(i); }
        const real_t *const* GetDiffusionCoeff22() const { return this->adterm->GetDiffusionCoeff22(); }
        const real_t *GetDiffusionCoeff22(const len_t i) const { return this->adterm->GetDiffusionCoeff22(i); }

        len_t GetNumberOfNonZerosPerRow() const;
        len_t GetNumberOfNonZerosPerRow_jac() const;
        PredeterminedParameter *GetPredetermined() { return this->predetermined; }
        /**
         * Returns 'true' if all terms of this equation are
         * 'PredeterminedParameter's.
         */
        bool IsPredetermined() const { return (predetermined != nullptr); }

        /**
         * This quantity can be evaluated if it consists of only either
         * a 'PredeterminedParameter' or one or more 'EvaluableEquationTerm'
         * objects.
         */
        bool IsEvaluable() const { return (adterm == nullptr && terms.size() == 0); }

        void RebuildTerms(const real_t, const real_t, UnknownQuantityHandler*);

        void SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix*);
        void SetJacobianBlockBC(const len_t uqtyId, const len_t derivId, Matrix*);
        void SetMatrixElements(Matrix*, real_t*);
        void SetVectorElements(real_t*, const real_t*);
    };
}

#endif/*_DREAM_FVM_EQUATION_HPP*/
