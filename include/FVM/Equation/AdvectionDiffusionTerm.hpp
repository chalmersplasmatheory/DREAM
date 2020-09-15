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

        enum AdvectionInterpolationCoefficient::adv_interpolation advectionInterpolationMethod_r  = AdvectionInterpolationCoefficient::AD_INTERP_CENTRED;
        enum AdvectionInterpolationCoefficient::adv_interpolation advectionInterpolationMethod_p1 = AdvectionInterpolationCoefficient::AD_INTERP_CENTRED;
        enum AdvectionInterpolationCoefficient::adv_interpolation advectionInterpolationMethod_p2 = AdvectionInterpolationCoefficient::AD_INTERP_CENTRED;
        
        real_t fluxLimiterDampingFactor = 1.0;

        // The following set of variables are used for dynamic damping of flux limiters
        bool withDynamicFluxLimiterDamping = true;
        real_t dampingWithIteration = 1.0, t_prev = -1.0, dt_prev = -1.0;
        len_t iteration=0;
    public:
        AdvectionDiffusionTerm(Grid *g)
            : AdvectionTerm(g, true), DiffusionTerm(g, true) {}

        void Add(AdvectionTerm*);
        void Add(DiffusionTerm*);

        virtual len_t GetNumberOfNonZerosPerRow() const override;
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override;

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override;
        virtual void ResetCoefficients() override;
        void RebuildInterpolationCoefficients(UnknownQuantityHandler*);


        virtual void SetJacobianBlock(
            const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t *x
        ) {
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
        
        // set the interpolation
        void SetAdvectionInterpolationMethod(
            AdvectionInterpolationCoefficient::adv_interpolation intp, 
            FVM::fluxGridType fgType, len_t id, real_t damping_factor 
        ){
            this->fluxLimiterDampingFactor = damping_factor;
            if(fgType == FLUXGRIDTYPE_RADIAL){
                this->advectionInterpolationMethod_r = intp; 
                this->deltar->SetUnknownId(id);
            } else if(fgType == FLUXGRIDTYPE_P1){
                this->advectionInterpolationMethod_p1 = intp;
                this->delta1->SetUnknownId(id);
            } else if(fgType == FLUXGRIDTYPE_P2){
                this->advectionInterpolationMethod_p2 = intp;
                this->delta2->SetUnknownId(id);
            } 
        }
        void SetAdvectionBoundaryConditions(
            fluxGridType fgType, AdvectionInterpolationCoefficient::adv_bc bc_lower, 
            AdvectionInterpolationCoefficient::adv_bc bc_upper
        ){
            if(fgType == FLUXGRIDTYPE_RADIAL)
                this->deltar->SetBoundaryConditions(bc_lower,bc_upper);
            else if(fgType == FLUXGRIDTYPE_P1)
                this->delta1->SetBoundaryConditions(bc_lower,bc_upper);
            else if(fgType == FLUXGRIDTYPE_P2)
                this->delta2->SetBoundaryConditions(bc_lower,bc_upper);
        }

    };
}

#endif/*_DREAM_FVM_ADVECTION_DIFFUSION_TERM_HPP*/
