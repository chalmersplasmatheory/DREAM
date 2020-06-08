#ifndef _DREAM_FVM_ADVECTION_DIFFUSION_TERM_HPP
#define _DREAM_FVM_ADVECTION_DIFFUSION_TERM_HPP

#include <vector>
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM::FVM {
    class AdvectionDiffusionTerm : public AdvectionTerm, public DiffusionTerm {
    public:
        enum advdiff_interpolation {
            AD_INTERP_CENTRED,
            AD_INTERP_BACKWARD,
            AD_INTERP_FORWARD
        };

    private:
        std::vector<AdvectionTerm*> advectionterms;
        std::vector<DiffusionTerm*> diffusionterms;

        enum advdiff_interpolation interpolationMethod = AD_INTERP_CENTRED;
        
    public:
        AdvectionDiffusionTerm(Grid *g, enum advdiff_interpolation intp=AD_INTERP_CENTRED)
            : AdvectionTerm(g, true), DiffusionTerm(g, true), interpolationMethod(intp) {}

        void Add(AdvectionTerm*);
        void Add(DiffusionTerm*);

        virtual len_t GetNumberOfNonZerosPerRow() const override;
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override;

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override;
        virtual void ResetCoefficients() override;
        void RebuildInterpolationCoefficients();
        void SetInterpolationCoefficientValues(const real_t);

        virtual void SetJacobianBlock(
            const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t *x
        ) {
            for (auto it = advectionterms.begin(); it != advectionterms.end(); it++)
                (*it)->SetJacobianBlock(uqtyId, derivId, jac, x);

            for (auto it = diffusionterms.begin(); it != diffusionterms.end(); it++)
                (*it)->SetJacobianBlock(uqtyId, derivId, jac, x);
        }
        virtual void SetMatrixElements(Matrix *mat, real_t *rhs) {
            if (this->advectionterms.size() > 0)
                this->AdvectionTerm::SetMatrixElements(mat, rhs);
            if (this->diffusionterms.size() > 0)
                this->DiffusionTerm::SetMatrixElements(mat, rhs);
        }
        virtual void SetVectorElements(real_t *vec, const real_t *x) {
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
