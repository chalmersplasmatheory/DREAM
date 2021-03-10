#ifndef _DREAM_FVM_ADVECTION_DIFFUSION_TERM_HPP
#define _DREAM_FVM_ADVECTION_DIFFUSION_TERM_HPP

#include <vector>
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM::FVM {
    class AdvectionDiffusionTerm : public AdvectionTerm, public DiffusionTerm {

    private:
        std::vector<AdvectionTerm*> advectionterms;
        std::vector<DiffusionTerm*> diffusionterms;

    public:
        AdvectionDiffusionTerm(Grid *g)
            : AdvectionTerm(g, true), DiffusionTerm(g, true) {}

        void Add(AdvectionTerm*);
        void Add(DiffusionTerm*);

        virtual len_t GetNumberOfNonZerosPerRow() const override;
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override;

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override;
        virtual void ResetCoefficients() override;

        const std::vector<AdvectionTerm*>& GetAdvectionTerms() const { return advectionterms; }
        const std::vector<DiffusionTerm*>& GetDiffusionTerms() const { return diffusionterms; }

        virtual const real_t *GetRadialJacobianInterpolationCoeffs() const override {
            if (this->advectionterms.size() > 0)
                return this->advectionterms[0]->GetRadialJacobianInterpolationCoeffs();
            else if (this->diffusionterms.size() > 0)
                return this->diffusionterms[0]->GetRadialJacobianInterpolationCoeffs();
            else
                return nullptr;
        }

        virtual void SetJacobianBlock(
            const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t *x
        ) override {
            this->interp_mode = AdvectionInterpolationCoefficient::AD_INTERP_MODE_JACOBIAN;
            // Set diagonal block (assuming constant coefficients)
            if (uqtyId == derivId) {
                if (this->advectionterms.size() > 0)
                    this->AdvectionTerm::SetMatrixElements(jac, nullptr);
                if (this->diffusionterms.size() > 0)
                    this->DiffusionTerm::SetMatrixElements(jac, nullptr);
            }

            // Handle any off-diagonal blocks and/or non-linear coefficients
            for (auto it = advectionterms.begin(); it != advectionterms.end(); it++)
                (*it)->SetJacobianBlock(uqtyId, derivId, jac, x);

            for (auto it = diffusionterms.begin(); it != diffusionterms.end(); it++)
                (*it)->SetJacobianBlock(uqtyId, derivId, jac, x);
        }
        virtual void SetMatrixElements(Matrix *mat, real_t *rhs) override {
            if (this->advectionterms.size() > 0)
                this->AdvectionTerm::SetMatrixElements(mat, rhs);
            if (this->diffusionterms.size() > 0)
                this->DiffusionTerm::SetMatrixElements(mat, rhs);
        }
        virtual void SetVectorElements(real_t *vec, const real_t *x) override {
            if (this->advectionterms.size() > 0)
                this->AdvectionTerm::SetVectorElements(vec, x);
            if (this->diffusionterms.size() > 0)
                this->DiffusionTerm::SetVectorElements(vec, x);
        }

        virtual void SaveCoefficientsSFile(const std::string&) override;
        virtual void SaveCoefficientsSFile(SFile*) override;
        
    };
}

#endif/*_DREAM_FVM_ADVECTION_DIFFUSION_TERM_HPP*/
