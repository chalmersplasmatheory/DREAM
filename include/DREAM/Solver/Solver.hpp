#ifndef _DREAM_SOLVER_HPP
#define _DREAM_SOLVER_HPP
/* Definition of the abstract base class 'Solver', which
 * defines the interface for all equation solvers in DREAM.
 */

#include <vector>
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class Solver {
    protected:
        FVM::UnknownQuantityHandler *unknowns;
        // List of equations associated with unknowns (owned by the 'EquationSystem')
        std::vector<UnknownQuantityEquation*> *unknown_equations;
        std::vector<len_t> nontrivial_unknowns;

    public:
        Solver(FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*);
        virtual ~Solver() {}

        void BuildJacobian(const real_t, const real_t, FVM::BlockMatrix*);
        void BuildMatrix(const real_t, const real_t, FVM::BlockMatrix*, real_t*);
        void BuildVector(const real_t, const real_t, real_t*);
        void RebuildTerms(const real_t, const real_t);

        virtual const real_t *GetSolution() const = 0;
        virtual void Initialize(const len_t, std::vector<len_t>&);

        // Apply the solver to stage 'n' of the non-linear
        // equation system. Returns 'true' if the stage is
        // considered to be converged after the solve.
        virtual void Solve(const real_t n, const real_t dt) = 0;
    };
}

#endif/*_DREAM_SOLVER_HPP*/
