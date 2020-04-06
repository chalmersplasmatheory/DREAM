#ifndef _DREAM_SOLVER_SNES_HPP
#define _DREAM_SOLVER_SNES_HPP

#include <vector>
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/MatrixInverter.hpp"

namespace DREAM {
    class SolverLinearlyImplicit : public Solver {
    private:
        // Matrix to store Jacobian in
        FVM::BlockMatrix* matrix = nullptr;
        // Vector to store source terms in
        Vec petsc_S;
        // Vector to store solution in
        Vec petsc_sol;

        FVM::MatrixInverter *inverter;

        real_t t, dt;

    protected:
        virtual void initialize_internal(const len_t, std::vector<len_t>&) override;

    public:
        SolverLinearlyImplicit(FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*);
        virtual ~SolverLinearlyImplicit();

        virtual void Initialize(const len_t, std::vector<len_t>&);

        real_t CurrentTime() const { return this->t; }
        real_t CurrentTimeStep() const { return this->dt; }
        FVM::BlockMatrix *GetMatrix() { return this->matrix; }

        virtual void SetInitialGuess(const real_t*) override;
        virtual void Solve(const real_t, const real_t) override;
    };
}

#endif/*_DREAM_SOLVER_SNES_HPP*/
