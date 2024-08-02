#ifndef _DREAM_CONVERGENCE_CHECKER_HPP
#define _DREAM_CONVERGENCE_CHECKER_HPP

#include <string>
#include <unordered_map>
#include <vector>
#include "FVM/NormEvaluator.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/DiagonalPreconditioner.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"

namespace DREAM {
    class ConvergenceChecker : public FVM::NormEvaluator {
    private:
        std::unordered_map<len_t, real_t> absTols;
        std::unordered_map<len_t, real_t> relTols;

		bool ownsPreconditioner = false;
		DiagonalPreconditioner *precond=nullptr;
		std::vector<UnknownQuantityEquation*> *unknown_eqns=nullptr;

        void DefineAbsoluteTolerances();
        real_t GetDefaultAbsTol(const std::string&);

        len_t dx_size = 0;
        real_t *dx_buffer=nullptr;

        len_t nNontrivials;
        real_t *x_2norm=nullptr;
        real_t *dx_2norm=nullptr;

		bool saveConvergenceInfo = false;
		std::unordered_map<len_t, std::vector<std::vector<real_t>>> convergence_x;
		std::unordered_map<len_t, std::vector<std::vector<real_t>>> convergence_dx;

		std::unordered_map<len_t, std::vector<bool>> residual_conv;
		std::unordered_map<len_t, std::vector<real_t>> residual_conv_maxerr;

    public:
        ConvergenceChecker(
            FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*,
			const std::vector<len_t>&,
			DiagonalPreconditioner *precond=nullptr, const real_t reltol=1e-6,
			bool saveConvergenceInfo=false
        );
        virtual ~ConvergenceChecker();

        void AllocateBuffer(const len_t);
        void DeallocateBuffer();

        const std::string &GetNonTrivialName(const len_t i)
        { return this->unknowns->GetUnknown(nontrivials[i])->GetName(); }

        bool IsConverged(const real_t*, const real_t*, const len_t, bool verbose=false);
        bool IsConverged(const real_t*, const real_t*, const real_t*, const len_t, bool verbose=false);

		bool IsResidualConverged(const len_t, const real_t, const real_t*, bool);

        const real_t *GetErrorNorms() { return this->dx_2norm; }
        const real_t GetErrorScale(const len_t);

		void SaveData(SFile*, const std::string&);

        void SetAbsoluteTolerance(const len_t, const real_t);
        void SetRelativeTolerance(const real_t);
        void SetRelativeTolerance(const len_t, const real_t);

		void SaveConvergenceInfo(const len_t, const len_t, const real_t, const real_t);
		void SetResidualConverged(const len_t, const len_t, const bool, const real_t);
    };
}

#endif/*_DREAM_CONVERGENCE_CHECKER_HPP*/
