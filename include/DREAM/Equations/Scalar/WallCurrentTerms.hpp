#ifndef _DREAM_EQUATION_SCALAR_WALL_CURRENT_TERMS_HPP
#define _DREAM_EQUATION_SCALAR_WALL_CURRENT_TERMS_HPP

#include "FVM/Equation/EvaluableEquationTerm.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

/**
 * Implementation of a class which represents the j_||/(B/Bmin) term in Ampere's law.
 */
namespace DREAM {
    class PoloidalFluxAtEdgeTerm : public FVM::EvaluableEquationTerm {
    public:
        PoloidalFluxAtEdgeTerm(FVM::Grid* g) : FVM::EvaluableEquationTerm(g){}
        
        virtual len_t GetNumberOfNonZerosPerRow() const override{ return 0; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override{ return 0; }

        virtual real_t* Evaluate(real_t*, const real_t*, const len_t, const len_t) override {return nullptr;}
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}
        virtual void SetJacobianBlock(const len_t /*derivId*/, const len_t /*uqtyId*/, FVM::Matrix*, const real_t*) override {}
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override {}
        virtual void SetVectorElements(real_t*, const real_t*) override {}

    };
}

namespace DREAM {
    class SOLMutualInductanceTerm : public FVM::EvaluableEquationTerm {
    public:
        SOLMutualInductanceTerm(FVM::Grid* g) : FVM::EvaluableEquationTerm(g){}

        virtual len_t GetNumberOfNonZerosPerRow() const override{ return 0; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override{ return 0; }

        virtual real_t* Evaluate(real_t*, const real_t*, const len_t, const len_t) override {return nullptr;}
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}
        virtual void SetJacobianBlock(const len_t /*derivId*/, const len_t /*uqtyId*/, FVM::Matrix*, const real_t*) override {}
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override {}
        virtual void SetVectorElements(real_t*, const real_t*) override {}

    };
}


namespace DREAM {
    class TotalPlasmaCurrentFromJTot : public FVM::EvaluableEquationTerm {
    public:
        TotalPlasmaCurrentFromJTot(FVM::Grid* g) : FVM::EvaluableEquationTerm(g){}

        virtual len_t GetNumberOfNonZerosPerRow() const override{ return 0; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override{ return 0; }

        virtual real_t* Evaluate(real_t*, const real_t*, const len_t, const len_t) override {return nullptr;}
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}
        virtual void SetJacobianBlock(const len_t /*derivId*/, const len_t /*uqtyId*/, FVM::Matrix*, const real_t*) override {}
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override {}
        virtual void SetVectorElements(real_t*, const real_t*) override {}

    };
}

#endif /*_DREAM_EQUATION_SCALAR_WALL_CURRENT_TERMS_HPP*/
