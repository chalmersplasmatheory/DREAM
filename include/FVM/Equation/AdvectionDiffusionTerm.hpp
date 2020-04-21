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
        AdvectionDiffusionTerm(Grid *g, enum advdiff_interpolation intp)
            : AdvectionTerm(g, true), DiffusionTerm(g, true), interpolationMethod(intp) {}

        void Add(AdvectionTerm*);
        void Add(DiffusionTerm*);

        virtual len_t GetNumberOfNonZerosPerRow() const override;
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override;

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override;
        void RebuildInterpolationCoefficients();
        void SetInterpolationCoefficientValues(const real_t);

        virtual void SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix *jac) {
            for (auto it = advectionterms.begin(); it != advectionterms.end(); it++)
                (*it)->SetJacobianBlock(uqtyId, derivId, jac);

            for (auto it = diffusionterms.begin(); it != diffusionterms.end(); it++)
                (*it)->SetJacobianBlock(uqtyId, derivId, jac);
        }
        virtual void SetMatrixElements(Matrix *mat, real_t *rhs) {
            this->AdvectionTerm::SetMatrixElements(mat, rhs);
            this->DiffusionTerm::SetMatrixElements(mat, rhs);
        }
        virtual void SetVectorElements(real_t *vec, const real_t *x) {
            this->AdvectionTerm::SetVectorElements(vec, x);
            this->DiffusionTerm::SetVectorElements(vec, x);
        }

        virtual void SaveCoefficientsSFile(const std::string&) override;
        virtual void SaveCoefficientsSFile(SFile*) override;
    };
}

#endif/*_DREAM_FVM_ADVECTION_DIFFUSION_TERM_HPP*/
