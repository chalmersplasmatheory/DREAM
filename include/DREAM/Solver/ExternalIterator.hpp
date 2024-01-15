#ifndef _DREAM_SOLVER_EXTERNAL_ITERATOR_HPP
#define _DREAM_SOLVER_EXTERNAL_ITERATOR_HPP

#include <vector>
#include "DREAM/ConvergenceChecker.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
	class ExternalIterator {
	protected:
		FVM::UnknownQuantityHandler *unknowns;

		std::vector<len_t> external_unknowns;
		std::vector<UnknownQuantityEquation*> *unknown_equations;
		len_t vecsize;
		real_t *vec=nullptr, *dvec=nullptr;

		bool printVerbose = false;

		ConvergenceChecker *convChecker=nullptr;
	
	public:
		ExternalIterator(
			FVM::UnknownQuantityHandler*,
			std::vector<UnknownQuantityEquation*>*,
			const bool verbose=false
		);
		~ExternalIterator();

		void Initialize(std::vector<len_t>&);
		bool Solve(const real_t, const real_t, const len_t);

		void SetConvergenceChecker(ConvergenceChecker*);
		void SetVerbose(bool v) { this->printVerbose = v; }
	};
}

#endif/*_DREAM_SOLVER_EXTERNAL_ITERATOR_HPP*/
