#ifndef _DREAM_CONVERGENCE_CHECKER_HPP
#define _DREAM_CONVERGENCE_CHECKER_HPP

#include <string>
#include <unordered_map>
#include <vector>
#include "FVM/NormEvaluator.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class ConvergenceChecker : public FVM::NormEvaluator {
    private:
        std::unordered_map<len_t, real_t> absTols;
        std::unordered_map<len_t, real_t> relTols;

        void DefineAbsoluteTolerances();
        real_t GetDefaultAbsTol(const std::string&);

        len_t dx_size = 0;
        real_t *dx_buffer=nullptr;

        len_t nNontrivials;
        real_t *x_2norm=nullptr;
        real_t *dx_2norm=nullptr;

    public:
        ConvergenceChecker(
            FVM::UnknownQuantityHandler*, const std::vector<len_t>&,
            const real_t reltol=1e-6
        );
        virtual ~ConvergenceChecker();

        void AllocateBuffer(const len_t);
        void DeallocateBuffer();

        const std::string &GetNonTrivialName(const len_t i)
        { return this->unknowns->GetUnknown(nontrivials[i])->GetName(); }

        bool IsConverged(const real_t*, const real_t*, bool verbose=false);
        bool IsConverged(const real_t*, const real_t*, const real_t*, bool verbose=false);

        const real_t *GetErrorNorms() { return this->dx_2norm; }
        const real_t GetErrorScale(const len_t);

        void SetAbsoluteTolerance(const len_t, const real_t);
        void SetRelativeTolerance(const real_t);
        void SetRelativeTolerance(const len_t, const real_t);
    };
}

#endif/*_DREAM_CONVERGENCE_CHECKER_HPP*/
