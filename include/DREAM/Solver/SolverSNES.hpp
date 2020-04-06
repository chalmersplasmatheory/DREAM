#ifndef _DREAM_SOLVER_SNES_HPP
#define _DREAM_SOLVER_SNES_HPP

#include <vector>
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/BlockMatrix.hpp"

namespace DREAM {
    class SolverSNES : public Solver {
    private:
        // Matrix to store Jacobian in
        FVM::BlockMatrix* jacobian = nullptr;
        // Vector to store function evaluations in
        Vec petsc_F;
        // Vector to store solution in
        Vec petsc_sol;

        SNES snes;

        real_t t, dt;

    public:
        SolverSNES(FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*);
        virtual ~SolverSNES();

        //virtual const real_t *GetSolution() const override { return F; }
        virtual void Initialize(const len_t, std::vector<len_t>&);

        real_t CurrentTime() const { return this->t; }
        real_t CurrentTimeStep() const { return this->dt; }
        FVM::BlockMatrix *GetJacobian() { return this->jacobian; }

        virtual void SetInitialGuess(const real_t*) override;
        virtual void Solve(const real_t, const real_t) override;
    };

    PetscErrorCode SNES_set_jacobian(SNES, Vec, Mat, Mat, void*);
    PetscErrorCode SNES_set_function(SNES, Vec, Vec, void*);
    void SNES_update_system(SolverSNES*);
    PetscErrorCode SNES_solution_obtained(SNES, PetscInt, PetscReal, void*);
}

#endif/*_DREAM_SOLVER_SNES_HPP*/
