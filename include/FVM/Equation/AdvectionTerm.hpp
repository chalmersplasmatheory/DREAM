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

        void AllocateCoefficients();
        void AllocateDifferentiationCoefficients();
        void AllocateInterpolationCoefficients();
        void DeallocateCoefficients();
        void DeallocateDifferentiationCoefficients();
        void DeallocateInterpolationCoefficients();        
        
        virtual void SetPartialAdvectionTerm(len_t /*derivId*/, len_t /*nMultiples*/){}
        void ResetJacobianColumn();
        std::vector<len_t> derivIds;
        std::vector<len_t> derivNMultiples;
        // Return maximum nMultiples for allocation of df
        len_t MaxNMultiple()
            {
            len_t nMultiples = 0;
            for(len_t it=0; it<derivIds.size(); it++)
                if (derivNMultiples[it]>nMultiples)
                    nMultiples = derivNMultiples[it];
            return nMultiples;
            }

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
/*
        const real_t *const* GetInterpolationCoeffR() const { return this->deltar; }
        const real_t *GetInterpolationCoeffR(const len_t i) const { return this->deltar[i]; }
        const real_t *const* GetInterpolationCoeff1() const { return this->delta1; }
        const real_t *GetInterpolationCoeff1(const len_t i) const { return this->delta1[i]; }
        const real_t *const* GetInterpolationCoeff2() const { return this->delta2; }
        const real_t *GetInterpolationCoeff2(const len_t i) const { return this->delta2[i]; }

        const real_t *GetInterpolationCoeffR(const len_t i) const { return this->deltar->GetCoefficient(i,0,0,1); }
        const real_t *GetInterpolationCoeff1(const len_t i) const { return this->delta1->GetCoefficient(0,i,0,1); }
        const real_t *GetInterpolationCoeff2(const len_t i) const { return this->delta2->GetCoefficient(0,0,i,1); }
*/
        const real_t GetInterpolationCoeff1(const len_t ir, const len_t i, const len_t j, const len_t n){return this->delta1->GetCoefficient(ir,i,j,n);}
        const real_t* GetInterpolationCoeff1(const len_t ir, const len_t i, const len_t j) const {return this->delta1->GetCoefficient(ir,i,j); }


        // TODO: FIX NNZ
        virtual len_t GetNumberOfNonZerosPerRow() const override 
            { return std::max(deltar->GetNNZPerRow(),std::max(delta1->GetNNZPerRow(),delta2->GetNNZPerRow())); }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override 
        { 
            len_t nnz = GetNumberOfNonZerosPerRow(); 
            for(len_t i = 0; i<derivIds.size(); i++)
                nnz += derivNMultiples[i];
            return nnz;
        }

        virtual void ResetCoefficients();
        virtual void ResetDifferentiationCoefficients();
        void SetCoefficients(
            real_t**, real_t**, real_t**, real_t**
        );

        // Accessors to advection coefficients
        real_t& Fr(const len_t ir, const len_t i1, const len_t i2)
        { return Fr(ir, i1, i2, this->fr); }
        real_t& Fr(const len_t ir, const len_t i1, const len_t i2, real_t **fr) {
            if (ir == nr) return fr[ir][i2*n1[ir-1] + i1];
            else return fr[ir][i2*n1[ir] + i1];
        }
        const real_t Fr(const len_t ir, const len_t i1, const len_t i2, const real_t *const* fr) const {
            if (ir == nr) return fr[ir][i2*n1[ir-1] + i1];
            else return fr[ir][i2*n1[ir] + i1];
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
            real_t*, const real_t*,
            const real_t *const*, const real_t *const*, const real_t *const*,const real_t *const*
        );

        // Adds derivId to list of unknown quantities that contributes to Jacobian of this advection term
        void AddUnknownForJacobian(FVM::UnknownQuantityHandler *u, len_t derivId){
            derivIds.push_back(derivId);
            derivNMultiples.push_back(u->GetUnknown(derivId)->NumberOfMultiples());
        }


        void SetInterpolationCoefficients(AdvectionInterpolationCoefficient*, AdvectionInterpolationCoefficient*, AdvectionInterpolationCoefficient*);

        virtual void SaveCoefficientsSFile(const std::string&);
        virtual void SaveCoefficientsSFile(SFile*);
    };
}

#endif/*_DREAM_FVM_ADVECTION_TERM_HPP*/
