#ifndef _DREAM_EQUATIONS_BRAAMSKARNEY_POTENTIAL_INTEGRAL_HPP
#define _DREAM_EQUATIONS_BRAAMSKARNEY_POTENTIAL_INTEGRAL_HPP

#include "FVM/Equation/EquationTerm.hpp"

namespace DREAM {

    class BraamsKarneyPotentialIntegral : public FVM::EquationTerm {
    public:
        enum class Type {
            UPS0, UPS1, UPS2, PI0, PI1
        };

        BraamsKarneyPotentialIntegral(FVM::Grid *grid, Type type);

        virtual len_t GetNumberOfNonZerosPerRow() const override {
            len_t nnz_per_row=0;
            for (len_t i = 0; i < grid->GetNr(); i++) {
                len_t nc = grid->GetMomentumGrid(i)->GetNCells();
                if (nnz_per_row < nc)
                    nnz_per_row = nc;
            }
            return nnz_per_row;
        }
        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}

        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

    private:
        Type type;
        std::map<std::tuple<len_t, len_t, len_t>, std::vector<real_t>> cache;

        template<typename F1, typename F2>
        void SetElements(F1 X, F2 ApplyX);

        std::tuple<std::vector<real_t>&, bool> getCachedIntegrand(len_t ir, len_t i, len_t j) {
            if (auto it {cache.find({ir, i, j})}; it != std::end(cache)) {
                return {it->second, true};
            } else {
                std::vector<real_t> &integrand = cache[{ir, i, j}];
                return {integrand, false};
            }
        }

    protected:
        real_t Integrand(real_t p, real_t xi, real_t pprime, real_t xiprime);

    };
}

#endif/*_DREAM_EQUATIONS_BRAAMSKARNEY_POTENTIAL_INTEGRAL_HPP*/
