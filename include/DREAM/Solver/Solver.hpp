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

        // Mapping from EquationSystem 'unknown_quantity_id' to index
        // in the block matrix representing the system
        std::map<len_t, len_t> unknownToMatrixMapping;

        // Number of rows in any (jacobian) matrix built by
        // this solver (not counting unknowns which should
        // not appear in the matrix)
        len_t matrix_size;

        virtual void initialize_internal(const len_t, std::vector<len_t>&) {}

    public:
        Solver(FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*);
        virtual ~Solver() {}

        void BuildJacobian(const real_t, const real_t, FVM::BlockMatrix*);
        void BuildMatrix(const real_t, const real_t, FVM::BlockMatrix*, real_t*);
        void BuildVector(const real_t, const real_t, real_t*);
        void RebuildTerms(const real_t, const real_t);

        //virtual const real_t *GetSolution() const = 0;
        virtual void Initialize(const len_t, std::vector<len_t>&);

        virtual void SetInitialGuess(const real_t*) = 0;
        virtual void Solve(const real_t t, const real_t dt) = 0;
    };
}

#endif/*_DREAM_SOLVER_HPP*/
