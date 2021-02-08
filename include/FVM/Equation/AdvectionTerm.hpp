#ifndef _DREAM_FVM_ADVECTION_TERM_HPP
#define _DREAM_FVM_ADVECTION_TERM_HPP

namespace DREAM::FVM { class AdvectionTerm; }

#include <softlib/SFile.h>
#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/AdvectionInterpolationCoefficient.hpp"
#include "FVM/Grid/Grid.hpp"

namespace DREAM::FVM {
    class AdvectionTerm : public EquationTerm {

    protected:

        // this is invoked when setting off-diagonal jacobian 
        // contributions on the radial flux grid 
        enum jacobian_interp_mode {
            NO_JACOBIAN         = 1, // used to SetVector and SetMatrixElements
            JACOBIAN_SET_LOWER  = 2, // sets offset contribution for radial flux  
            JACOBIAN_SET_CENTER = 3, // sets diagonal jacobian contributions
            JACOBIAN_SET_UPPER  = 4  // sets offset contribution for radial flux
        };
        // interpolation coefficients to radial flux grid
        real_t *deltaRadialFlux=nullptr; 

        real_t 
            **fr = nullptr, 
            **f1 = nullptr, 
            **f2 = nullptr;
        real_t
            **dfr = nullptr,
            **df1 = nullptr,
            **df2 = nullptr;
        real_t **f1pSqAtZero   = nullptr;
        real_t **df1pSqAtZero  = nullptr;
        real_t *JacobianColumn = nullptr;

        bool coefficientsShared = false;

        // Interpolation coefficients
        AdvectionInterpolationCoefficient *deltar=nullptr, *delta1=nullptr, *delta2=nullptr;
        bool interpolationCoefficientsShared = false;

        enum AdvectionInterpolationCoefficient::adv_interpolation advectionInterpolationMethod_r  = AdvectionInterpolationCoefficient::AD_INTERP_CENTRED;
        enum AdvectionInterpolationCoefficient::adv_interpolation advectionInterpolationMethod_p1 = AdvectionInterpolationCoefficient::AD_INTERP_CENTRED;
        enum AdvectionInterpolationCoefficient::adv_interpolation advectionInterpolationMethod_p2 = AdvectionInterpolationCoefficient::AD_INTERP_CENTRED;
        
        real_t fluxLimiterDampingFactor = 1.0;

        // The following set of variables are used for dynamic damping of flux limiters
        bool withDynamicFluxLimiterDamping = true;
        real_t dampingWithIteration = 1.0, t_prev = -1.0, dt_prev = -1.0;
        len_t iteration=0;

        // Memory allocation for various coefficients
        void AllocateCoefficients();
        void AllocateDifferentiationCoefficients();
        void AllocateInterpolationCoefficients();
        void DeallocateCoefficients();
        void DeallocateDifferentiationCoefficients();
        void DeallocateInterpolationCoefficients();        
        
        void SetPartialJacobianContribution(int_t, jacobian_interp_mode, len_t, Matrix*, const real_t*);
        void ResetJacobianColumn();

        AdvectionInterpolationCoefficient::adv_interp_mode interp_mode
            = AdvectionInterpolationCoefficient::AD_INTERP_MODE_FULL;
    public:
        AdvectionTerm(Grid*, bool allocateCoeffs=false);
        ~AdvectionTerm();

        const real_t *const* GetAdvectionCoeffR() const { return this->fr; }
        const real_t *GetAdvectionCoeffR(const len_t i) const { return this->fr[i]; }
        const real_t *const* GetAdvectionCoeff1() const { return this->f1; }
        const real_t *GetAdvectionCoeff1(const len_t i) const { return this->f1[i]; }
        const real_t *const* GetAdvectionCoeff2() const { return this->f2; }
        const real_t *GetAdvectionCoeff2(const len_t i) const { return this->f2[i]; }

        const real_t *const* GetAdvectionDiffCoeffR() const { return this->dfr; }
        const real_t *GetAdvectionDiffCoeffR(const len_t i) const { return this->dfr[i]; }
        const real_t *const* GetAdvectionDiffCoeff1() const { return this->df1; }
        const real_t *GetAdvectionDiffCoeff1(const len_t i) const { return this->df1[i]; }
        const real_t *const* GetAdvectionDiffCoeff2() const { return this->df2; }
        const real_t *GetAdvectionDiffCoeff2(const len_t i) const { return this->df2[i]; }
        
        virtual const real_t *GetRadialJacobianInterpolationCoeffs() const { return deltaRadialFlux; }

        // Returns nnz per row (assuming that this AdvectionTerm contains non-zero
        // elements in all three components)
        virtual len_t GetNumberOfNonZerosPerRow() const override {
            return 1 +
                deltar->GetOffDiagonalNNZPerRow() +
                delta1->GetOffDiagonalNNZPerRow() +
                delta2->GetOffDiagonalNNZPerRow();
        }

        virtual void ResetCoefficients();
        virtual void ResetDifferentiationCoefficients();
        void SetCoefficients(
            real_t**, real_t**, real_t**, real_t**, real_t*
        );

        // Accessors to advection coefficients
        real_t& Fr(const len_t ir, const len_t i1, const len_t i2)
        { return Fr(ir, i1, i2, this->fr); }
        real_t& Fr(const len_t ir, const len_t i1, const len_t i2, real_t **fr) {
            len_t np1 = (ir==nr) ? n1[ir-1] : n1[ir];
            return fr[ir][i2*np1 + i1];
        }
        const real_t Fr(const len_t ir, const len_t i1, const len_t i2, const real_t *const* fr) const {
            len_t np1 = (ir==nr) ? n1[ir-1] : n1[ir];
            return fr[ir][i2*np1 + i1];
        }

        real_t& F1(const len_t ir, const len_t i1, const len_t i2)
        { return F1(ir, i1, i2, this->f1); }
        real_t& F1(const len_t ir, const len_t i1, const len_t i2, real_t **f1)
        { return f1[ir][i2*(n1[ir]+1) + i1]; }
        const real_t F1(const len_t ir, const len_t i1, const len_t i2, const real_t *const* f1) const
        { return f1[ir][i2*(n1[ir]+1) + i1]; }

        real_t& F2(const len_t ir, const len_t i1, const len_t i2)
        { return F2(ir, i1, i2, this->f2); }
        real_t& F2(const len_t ir, const len_t i1, const len_t i2, real_t **f2)
        { return f2[ir][i2*n1[ir] + i1]; }
        const real_t F2(const len_t ir, const len_t i1, const len_t i2, const real_t *const* f2) const
        { return f2[ir][i2*n1[ir] + i1]; }

        real_t& F1PSqAtZero(const len_t ir, const len_t i2)
        { return F1PSqAtZero(ir, i2, this->f1pSqAtZero); }
        real_t& F1PSqAtZero(const len_t ir, const len_t i2, real_t **f1pSqAtZero)
        { return f1pSqAtZero[ir][i2]; }
        const real_t F1PSqAtZero(const len_t ir, const len_t i2, const real_t *const* f1pSqAtZero) const
        { return f1pSqAtZero[ir][i2]; }


        // Accessors to differentiation coefficients
        real_t& dFr(const len_t ir, const len_t i1, const len_t i2, const len_t nMultiple) {
            if (ir == nr) return dfr[ir+(nr+1)*nMultiple][i2*n1[ir-1] + i1];
            else return dfr[ir+(nr+1)*nMultiple][i2*n1[ir] + i1];
        }
        real_t& dF1(const len_t ir, const len_t i1, const len_t i2, const len_t nMultiple)
        { return df1[ir+nr*nMultiple][i2*(n1[ir]+1) + i1]; }
        real_t& dF2(const len_t ir, const len_t i1, const len_t i2, const len_t nMultiple)
        { return df2[ir+nr*nMultiple][i2*n1[ir] + i1]; }
        real_t& dF1PSqAtZero(const len_t ir, const len_t i2, const len_t nMultiple)
        { return df1pSqAtZero[ir+nr*nMultiple][i2]; }

        virtual bool GridRebuilt() override;
        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual void SetVectorElements(
            real_t*, const real_t*, const real_t *const*, 
            const real_t *const*, const real_t *const*,const real_t *const*, jacobian_interp_mode set=NO_JACOBIAN
        );

        void SetInterpolationCoefficients(AdvectionInterpolationCoefficient*, AdvectionInterpolationCoefficient*, AdvectionInterpolationCoefficient*);

        const real_t GetInterpolationCoeffR(const len_t ir, const len_t i, const len_t j, const len_t n) const { return this->deltar->GetCoefficient(ir, i, j, n); }
        const real_t *GetInterpolationCoeffR(const len_t ir, const len_t i, const len_t j) const { return this->deltar->GetCoefficient(ir, i, j); }
        const real_t GetInterpolationCoeff1(const len_t ir, const len_t i, const len_t j, const len_t n) const { return this->delta1->GetCoefficient(ir, i, j, n); }
        const real_t *GetInterpolationCoeff1(const len_t ir, const len_t i, const len_t j) const { return this->delta1->GetCoefficient(ir, i, j); }
        const real_t GetInterpolationCoeff2(const len_t ir, const len_t i, const len_t j, const len_t n) const { return this->delta2->GetCoefficient(ir, i, j, n); }
        const real_t *GetInterpolationCoeff2(const len_t ir, const len_t i, const len_t j) const { return this->delta2->GetCoefficient(ir, i, j); }

        AdvectionInterpolationCoefficient *GetInterpolationCoeffR() { return this->deltar; }
        AdvectionInterpolationCoefficient *GetInterpolationCoeff1() { return this->delta1; }
        AdvectionInterpolationCoefficient *GetInterpolationCoeff2() { return this->delta2; }

        void RebuildFluxLimiterDamping(const real_t, const real_t);
        void RebuildInterpolationCoefficients(UnknownQuantityHandler*, real_t**, real_t**, real_t**);

        virtual void SaveCoefficientsSFile(const std::string&);
        virtual void SaveCoefficientsSFile(SFile*);

        virtual void SetPartialAdvectionTerm(len_t /*derivId*/, len_t /*nMultiples*/){}

        // set the interpolation
        void SetAdvectionInterpolationMethod(
            AdvectionInterpolationCoefficient::adv_interpolation intp,
            OptionConstants::adv_jacobian_mode jac_mode, 
            FVM::fluxGridType fgType, len_t id, real_t damping_factor 
        ){
            this->fluxLimiterDampingFactor = damping_factor;
            if(fgType == FLUXGRIDTYPE_RADIAL){
                this->advectionInterpolationMethod_r = intp; 
                this->deltar->SetUnknownId(id);
                this->deltar->SetJacobianMode(jac_mode);
            } else if(fgType == FLUXGRIDTYPE_P1){
                this->advectionInterpolationMethod_p1 = intp;
                this->delta1->SetUnknownId(id);
                this->delta1->SetJacobianMode(jac_mode);
            } else if(fgType == FLUXGRIDTYPE_P2){
                this->advectionInterpolationMethod_p2 = intp;
                this->delta2->SetUnknownId(id);
                this->delta2->SetJacobianMode(jac_mode);
            } 
        }
        // set same interpolation method on all components 
        void SetAdvectionInterpolationMethod(
            AdvectionInterpolationCoefficient::adv_interpolation intp,
            OptionConstants::adv_jacobian_mode jac_mode, 
            len_t id, real_t damping_factor=1.0 
        ){
            this->fluxLimiterDampingFactor = damping_factor;
            this->advectionInterpolationMethod_r = intp; 
            this->deltar->SetUnknownId(id);
            this->deltar->SetJacobianMode(jac_mode);
            this->advectionInterpolationMethod_p1 = intp;
            this->delta1->SetUnknownId(id);
            this->delta1->SetJacobianMode(jac_mode);
            this->advectionInterpolationMethod_p2 = intp;
            this->delta2->SetUnknownId(id);
            this->delta2->SetJacobianMode(jac_mode);
        }

        // set boundary conditions
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
        // set same boundary conditions on all components
        void SetAdvectionBoundaryConditions(
            AdvectionInterpolationCoefficient::adv_bc bc_lower, 
            AdvectionInterpolationCoefficient::adv_bc bc_upper
        ){
            this->deltar->SetBoundaryConditions(bc_lower,bc_upper);
            this->delta1->SetBoundaryConditions(bc_lower,bc_upper);
            this->delta2->SetBoundaryConditions(bc_lower,bc_upper);
        }
    };
}

#endif/*_DREAM_FVM_ADVECTION_TERM_HPP*/
