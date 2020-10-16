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

        virtual void SetPartialDiffusionTerm(len_t /*derivId*/, len_t /*nMultiples*/){}
        void ResetJacobianColumn();
        std::vector<len_t> derivIds;
        std::vector<len_t> derivNMultiples;
        
        // Return maximum nMultiples for allocation of dd
        len_t MaxNMultiple() {
            len_t nMultiples = 0;
            for(len_t it=0; it<derivIds.size(); it++)
                if (derivNMultiples[it]>nMultiples)
                    nMultiples = derivNMultiples[it];
            return nMultiples;
        }


    public:
        DiffusionTerm(Grid*, bool allocCoefficients=false);
        virtual ~DiffusionTerm();

        void AllocateCoefficients();
        void AllocateDifferentiationCoefficients();
        void DeallocateCoefficients();
        void DeallocateDifferentiationCoefficients();
        void SetCoefficients(
            real_t**, real_t**, real_t**, real_t**, real_t**
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

        virtual len_t GetNumberOfNonZerosPerRow() const override {
            len_t nnz = 1;

            len_t np1 = this->grid->GetMomentumGrid(0)->GetNp1();
            len_t np2 = this->grid->GetMomentumGrid(0)->GetNp2();

            if (this->grid->GetNr() > 1) nnz += 2;      // Drr
            // XXX here we assume that all momentum grids are the same
            if (np1 > 1) nnz += 2;      // Dpp
            if (np2 > 1) nnz += 2;      // Dxx
            if (np1 > 1 && np2 > 1) nnz += 4;   // Dpx & Dxp

            return nnz;
        }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { 
                len_t nnz = GetNumberOfNonZerosPerRow(); 
                for(len_t i = 0; i<derivIds.size(); i++)
                    nnz += derivNMultiples[i];
                return nnz;
            }
            
        // Accessors to diffusion coefficients
        real_t& Drr(const len_t ir, const len_t i1, const len_t i2)
        { return Drr(ir, i1, i2, this->drr); }
        // XXX here we explicitly assume that the momentum
        // grids are the same at all radii
        real_t& Drr(const len_t ir, const len_t i1, const len_t i2, real_t **drr) {
            if (ir == nr) return drr[ir][i2*n1[ir-1] + i1];
            else return drr[ir][i2*n1[ir] + i1];
        }
        const real_t Drr(const len_t ir, const len_t i1, const len_t i2, const real_t *const* drr) const {
            if (ir == nr) return drr[ir][i2*n1[ir-1] + i1];
            else return drr[ir][i2*n1[ir] + i1];
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
            if (ir == nr)
                // XXX here we explicitly assume that the momentum
                // grids are the same at all radii
                return ddrr[ir+nMultiple*(nr+1)][i2*n1[ir-1] + i1];
            else
                return ddrr[ir+nMultiple*(nr+1)][i2*n1[ir] + i1];
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
            const real_t *const*, const real_t *const*
        );

        // Adds derivId to list of unknown quantities that contributes to Jacobian of this diffusion term
        void AddUnknownForJacobian(FVM::UnknownQuantityHandler *u, len_t derivId){
            derivIds.push_back(derivId);
            derivNMultiples.push_back(u->GetUnknown(derivId)->NumberOfMultiples());
        }


        virtual void SaveCoefficientsSFile(const std::string&);
        virtual void SaveCoefficientsSFile(SFile*);
    };
}


#endif/*_DREAM_FVM_DIFFUSION_TERM_HPP*/
