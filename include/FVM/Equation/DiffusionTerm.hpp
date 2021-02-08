namespace DREAM::FVM { class DiffusionTerm; }

#include <softlib/SFile.h>
#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"

#ifndef _DREAM_FVM_DIFFUSION_TERM_HPP
#define _DREAM_FVM_DIFFUSION_TERM_HPP

namespace DREAM::FVM {
    class DiffusionTerm : public EquationTerm {
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
            **drr=nullptr,
            **d11=nullptr, **d12=nullptr,
            **d21=nullptr, **d22=nullptr;
        real_t
            **ddrr=nullptr,
            **dd11=nullptr, **dd12=nullptr,
            **dd21=nullptr, **dd22=nullptr;
        real_t *JacobianColumn = nullptr;


        bool coefficientsShared = false;

        void SetPartialJacobianContribution(int_t, jacobian_interp_mode, len_t, Matrix*, const real_t*);
        void ResetJacobianColumn();

    public:
        DiffusionTerm(Grid*, bool allocCoefficients=false);
        virtual ~DiffusionTerm();

        void AllocateCoefficients();
        void AllocateDifferentiationCoefficients();
        void DeallocateCoefficients();
        void DeallocateDifferentiationCoefficients();
        void SetCoefficients(
            real_t**, real_t**, real_t**, real_t**, real_t**, real_t*
        );
        virtual void ResetCoefficients();
        virtual void ResetDifferentiationCoefficients();

        const real_t *const* GetDiffusionCoeffRR() const { return this->drr; }
        const real_t *GetDiffusionCoeffRR(const len_t i) const { return this->drr[i]; }
        const real_t *const* GetDiffusionCoeff11() const { return this->d11; }
        const real_t *GetDiffusionCoeff11(const len_t i) const { return this->d11[i]; }
        const real_t *const* GetDiffusionCoeff12() const { return this->d12; }
        const real_t *GetDiffusionCoeff12(const len_t i) const { return this->d12[i]; }
        const real_t *const* GetDiffusionCoeff21() const { return this->d21; }
        const real_t *GetDiffusionCoeff21(const len_t i) const { return this->d21[i]; }
        const real_t *const* GetDiffusionCoeff22() const { return this->d22; }
        const real_t *GetDiffusionCoeff22(const len_t i) const { return this->d22[i]; }

        const real_t *const* GetDiffusionDiffCoeffRR() const { return this->ddrr; }
        const real_t *GetDiffusionDiffCoeffRR(const len_t i) const { return this->ddrr[i]; }
        const real_t *const* GetDiffusionDiffCoeff11() const { return this->dd11; }
        const real_t *GetDiffusionDiffCoeff11(const len_t i) const { return this->dd11[i]; }
        const real_t *const* GetDiffusionDiffCoeff12() const { return this->dd12; }
        const real_t *GetDiffusionDiffCoeff12(const len_t i) const { return this->dd12[i]; }
        const real_t *const* GetDiffusionDiffCoeff21() const { return this->dd21; }
        const real_t *GetDiffusionDiffCoeff21(const len_t i) const { return this->dd21[i]; }
        const real_t *const* GetDiffusionDiffCoeff22() const { return this->dd22; }
        const real_t *GetDiffusionDiffCoeff22(const len_t i) const { return this->dd22[i]; }

        virtual len_t GetNumberOfNonZerosPerRow() const override {
            len_t nnz = 1;

            len_t np1 = this->grid->GetMomentumGrid(0)->GetNp1();
            len_t np2 = this->grid->GetMomentumGrid(0)->GetNp2();

            if (this->grid->GetNr() > 1) nnz += 2;      // Drr
            // XXX here we assume that all momentum grids are the same
            if (np1 > 1) nnz += 2;      // Dpp
            if (np2 > 1) nnz += 2;      // Dxx
            if (np1 > 1 && np2 > 1) nnz += 2*8;   // Dpx & Dxp

            return nnz;
        }

        // Accessors to diffusion coefficients
        real_t& Drr(const len_t ir, const len_t i1, const len_t i2)
        { return Drr(ir, i1, i2, this->drr); }
        // XXX here we explicitly assume that the momentum
        // grids are the same at all radii
        real_t& Drr(const len_t ir, const len_t i1, const len_t i2, real_t **drr) {
            len_t np1 = (ir==nr) ? n1[ir-1] : n1[ir];
            return drr[ir][i2*np1 + i1];
        }
        const real_t Drr(const len_t ir, const len_t i1, const len_t i2, const real_t *const* drr) const {
            len_t np1 = (ir==nr) ? n1[ir-1] : n1[ir];
            return drr[ir][i2*np1 + i1];
        }

        real_t& D11(const len_t ir, const len_t i1_f, const len_t i2)
        { return D11(ir, i1_f, i2, this->d11); }
        real_t& D11(const len_t ir, const len_t i1_f, const len_t i2, real_t **d11)
        { return d11[ir][i2*(n1[ir]+1) + i1_f]; }
        const real_t D11(const len_t ir, const len_t i1_f, const len_t i2, const real_t *const* d11) const
        { return d11[ir][i2*(n1[ir]+1) + i1_f]; }

        real_t& D12(const len_t ir, const len_t i1_f, const len_t i2)
        { return D12(ir, i1_f, i2, this->d12); }
        real_t& D12(const len_t ir, const len_t i1_f, const len_t i2, real_t **d12)
        { return d12[ir][i2*(n1[ir]+1) + i1_f]; }
        const real_t D12(const len_t ir, const len_t i1_f, const len_t i2, const real_t *const* d12) const
        { return d12[ir][i2*(n1[ir]+1) + i1_f]; }

        real_t& D21(const len_t ir, const len_t i1, const len_t i2_f)
        { return D21(ir, i1, i2_f, this->d21); }
        real_t& D21(const len_t ir, const len_t i1, const len_t i2_f, real_t **d21)
        { return d21[ir][i2_f*n1[ir] + i1]; }
        const real_t D21(const len_t ir, const len_t i1, const len_t i2_f, const real_t *const* d21) const
        { return d21[ir][i2_f*n1[ir] + i1]; }

        real_t& D22(const len_t ir, const len_t i1, const len_t i2_f)
        { return D22(ir, i1, i2_f, this->d22); }
        real_t& D22(const len_t ir, const len_t i1, const len_t i2_f, real_t **d22)
        { return d22[ir][i2_f*n1[ir] + i1]; }
        const real_t& D22(const len_t ir, const len_t i1, const len_t i2_f, const real_t *const* d22) const
        { return d22[ir][i2_f*n1[ir] + i1]; }

        real_t& dDrr(const len_t ir, const len_t i1, const len_t i2, const len_t nMultiple=0) {
            // XXX here we explicitly assume that the momentum
            // grids are the same at all radii
            len_t np1 = (ir==nr) ? n1[ir-1] : n1[ir];
            return ddrr[ir+nMultiple*(nr+1)][i2*np1 + i1];
        }
        real_t& dD11(const len_t ir, const len_t i1_f, const len_t i2, const len_t nMultiple=0)
        { return dd11[ir+nMultiple*nr][i2*(n1[ir]+1) + i1_f]; }
        real_t& dD12(const len_t ir, const len_t i1_f, const len_t i2, const len_t nMultiple=0)
        { return dd12[ir+nMultiple*nr][i2*(n1[ir]+1) + i1_f]; }
        real_t& dD21(const len_t ir, const len_t i1, const len_t i2_f, const len_t nMultiple=0)
        { return dd21[ir+nMultiple*nr][i2_f*n1[ir] + i1]; }
        real_t& dD22(const len_t ir, const len_t i1, const len_t i2_f, const len_t nMultiple=0)
        { return dd22[ir+nMultiple*nr][i2_f*n1[ir] + i1]; }

        virtual bool GridRebuilt() override;
        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual void SetVectorElements(
            real_t*, const real_t*,
            const real_t *const*, const real_t *const*, const real_t *const*,
            const real_t *const*, const real_t *const*, jacobian_interp_mode set=NO_JACOBIAN
        );

        // In this method the calculation of ddrr, dd11, ... should be implemented
        // in the derived classes, which represent derivatives of diffusion coefficients 
        // with respect to the unknown derivId. __We assume that the diffusion 
        // coefficient is local__, so that dxx in grid point ir only depends on the 
        // unknown in point ir, therefore ddxx will be of the same size as dxx (times nMultiples
        // when there are multiple such). For ddrr, the derivatives are taken with respect
        // to the unknown evaluated on the radial flux grid, which is interpolated to via
        // the 2-point stencil given by deltaRadialFlux.
        virtual void SetPartialDiffusionTerm(len_t /*derivId*/, len_t /*nMultiples*/){}

        virtual void SaveCoefficientsSFile(const std::string&);
        virtual void SaveCoefficientsSFile(SFile*);
    };
}


#endif/*_DREAM_FVM_DIFFUSION_TERM_HPP*/
