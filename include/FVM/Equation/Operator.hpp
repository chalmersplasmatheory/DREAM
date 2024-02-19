#ifndef _DREAM_FVM_EQUATION_HPP
#define _DREAM_FVM_EQUATION_HPP

#include <vector>
#include <unordered_map>
#include "FVM/Equation/AdvectionDiffusionTerm.hpp"
#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/EvaluableEquationTerm.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Grid/Grid.hpp"

namespace DREAM::FVM {
    class OperatorException : public FVMException {
    public:
        template<typename ... Args>
        OperatorException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("Operator");
        }
    };

    class Operator {
    private:
        std::vector<BC::BoundaryCondition*> boundaryConditions;
        std::vector<EquationTerm*> terms;
        std::vector<EvaluableEquationTerm*> eval_terms;
        PredeterminedParameter *predetermined = nullptr;
        AdvectionDiffusionTerm *adterm = nullptr;
        Grid *grid;
		bool hasTransientTerm = false;

        // List of pointers to terms which can be identified with a
        // numeric value and returned separately...
        std::unordered_map<int_t, EquationTerm*> identifiableTerms;

        real_t *vectorElementsSingleTerm=nullptr;

        void MakeIdentifiable(int_t, EquationTerm*);

    public:
        Operator(Grid*);

        ~Operator();

        void AddTerm(AdvectionTerm*, bool addAsEquationTerm = false);
        void AddTerm(DiffusionTerm*, bool addAsEquationTerm = false);
        void AddTerm(PredeterminedParameter*);
        void AddTerm(EvaluableEquationTerm*);
        void AddTerm(EquationTerm*);
		void AddTerm(TransientTerm*);
        void AddBoundaryCondition(BC::BoundaryCondition*);

        // Verifies that the operator is consistent
        void CheckConsistency();

        void Evaluate(real_t*, const real_t*);
        void EvaluableTransform(real_t*);

		bool HasTransientTerm() const { return this->hasTransientTerm; }

        EquationTerm *GetTermByID(const int_t id) { return this->identifiableTerms[id]; }

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

        bool IsEmpty() const;

/*
        const real_t *GetInterpolationCoeffR(const len_t i) const { return this->adterm->GetInterpolationCoeffR(i); }
        const real_t *GetInterpolationCoeff1(const len_t i) const { return this->adterm->GetInterpolationCoeff1(i); }
        const real_t *GetInterpolationCoeff2(const len_t i) const { return this->adterm->GetInterpolationCoeff2(i); }
*/
        AdvectionInterpolationCoefficient *GetInterpolationCoeffR() const { return this->adterm->GetInterpolationCoeffR(); }
        const real_t GetInterpolationCoeffR(const len_t ir, const len_t i, const len_t j, const len_t n) const { return this->adterm->GetInterpolationCoeffR(ir,i,j,n); }
        const real_t* GetInterpolationCoeffR(const len_t ir, const len_t i, const len_t j) const { return this->adterm->GetInterpolationCoeffR(ir,i,j); }
        AdvectionInterpolationCoefficient *GetInterpolationCoeff1() const { return this->adterm->GetInterpolationCoeff1(); }
        const real_t GetInterpolationCoeff1(const len_t ir, const len_t i, const len_t j, const len_t n) const { return this->adterm->GetInterpolationCoeff1(ir,i,j,n); }
        const real_t* GetInterpolationCoeff1(const len_t ir, const len_t i, const len_t j) const { return this->adterm->GetInterpolationCoeff1(ir,i,j); }
        AdvectionInterpolationCoefficient *GetInterpolationCoeff2() const { return this->adterm->GetInterpolationCoeff2(); }
        const real_t GetInterpolationCoeff2(const len_t ir, const len_t i, const len_t j, const len_t n) const { return this->adterm->GetInterpolationCoeff2(ir,i,j,n); }
        const real_t* GetInterpolationCoeff2(const len_t ir, const len_t i, const len_t j) const { return this->adterm->GetInterpolationCoeff2(ir,i,j); }

        void SetAdvectionInterpolationMethod(AdvectionInterpolationCoefficient::adv_interpolation intp, OptionConstants::adv_jacobian_mode jac, FVM::fluxGridType fgType, len_t id=0, real_t damping=1.0) 
            { this->adterm->SetAdvectionInterpolationMethod(intp,jac,fgType,id,damping); }
        void SetAdvectionBoundaryConditions(fluxGridType fluxGridType, AdvectionInterpolationCoefficient::adv_bc bc_lower, AdvectionInterpolationCoefficient::adv_bc bc_upper)
            { this->adterm->SetAdvectionBoundaryConditions(fluxGridType,bc_lower,bc_upper); }
        void SetAdvectionBoundaryConditions(AdvectionInterpolationCoefficient::adv_bc bc_lower, AdvectionInterpolationCoefficient::adv_bc bc_upper)
            { this->adterm->SetAdvectionBoundaryConditions(bc_lower,bc_upper); }

        len_t GetNumberOfNonZerosPerRow() const;
        len_t GetNumberOfNonZerosPerRow_jac() const;
        PredeterminedParameter *GetPredetermined() { return this->predetermined; }
        AdvectionDiffusionTerm *GetAdvectionDiffusion() const { return this->adterm; }
        /**
         * Returns 'true' if all terms of this operator are
         * 'PredeterminedParameter's.
         */
        bool IsPredetermined() const { return (predetermined != nullptr); }
        bool IsEvaluable() const;

        void RebuildTerms(const real_t, const real_t, UnknownQuantityHandler*);

        bool SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix*, const real_t*, bool printTerms=false);
        bool SetJacobianBlockBC(const len_t uqtyId, const len_t derivId, Matrix*, const real_t*, bool printTerms=false);
        void SetMatrixElements(Matrix*, real_t*);
        void SetVectorElements(real_t*, const real_t*);

        /**
         * Returns the contribution from the SetVectorElements-function from an individual term
         * To be used for accessing contributions from individual terms 
         * so they can be saved to the output as "other quantities"
         * TODO: come up with a better way to handle id's for individual terms, 
         * now one has to check in which order the terms are added to the operator!
         */ 
        const real_t* GetVectorElementsSingleEquationTerm(len_t, const real_t*) const;
    };
}

#endif/*_DREAM_FVM_EQUATION_HPP*/
