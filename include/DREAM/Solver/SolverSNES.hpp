#ifndef _DREAM_SOLVER_SNES_HPP
#define _DREAM_SOLVER_SNES_HPP

#include <vector>
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/BlockMatrix.hpp"

// Define potentially missing PETSc symbols
#if PETSC_VERSION_MAJOR < 3 || PETSC_VERSION_MINOR < 12
#   define SNES_DIVERGED_TR_DELTA SNES_CONVERGED_TR_DELTA
#endif

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

        // 2-norms of the solution vector and its residual
        real_t *x_2norm, *dx_2norm;

        // Maximum number of iterations allowed
        PetscInt maxiter;
        // Relative tolerance used for convergence test
        real_t reltol;

        // If true, this solver object will output detailed
        // information about the solver progress to stdout
        bool verbose;


        // Number of the iteration at which the equation
        // system was last rebuilt
        int_t last_rebuild = -1;
        
        void _EvaluateF(const real_t*, real_t*, FVM::BlockMatrix*);
        void _EvaluateJacobianNumerically(FVM::BlockMatrix*);

    protected:
        virtual void initialize_internal(const len_t, std::vector<len_t>&);

    public:
        SolverSNES(
            FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*,
            const PetscInt maxiter=100, const real_t reltol=1e-5,
            bool verbose=true
        );
        virtual ~SolverSNES();

        //virtual const real_t *GetSolution() const override { return F; }

        real_t CurrentTime() const { return this->t; }
        real_t CurrentTimeStep() const { return this->dt; }
        FVM::BlockMatrix *GetJacobian() { return this->jacobian; }
        int_t GetLastRebuild() const { return this->last_rebuild; }
        void SetLastRebuild(int_t i) { this->last_rebuild = i; }

        len_t GetNUnknowns() const { return this->nontrivial_unknowns.size(); }
        Vec GetPETScSolution() { return this->petsc_sol; }
        Vec GetPETScFunction() { return this->petsc_F; }

        // "Internal methods"
        const std::string& GetNonTrivialName(const len_t);
        real_t *Get_x_2norm() { return this->x_2norm; }
        real_t *Get_dx_2norm() { return this->dx_2norm; }

        PetscInt MaxIter() const { return this->maxiter; }
        real_t RelTol() const { return this->reltol; }
        bool Verbose() const { return this->verbose; }

        virtual void SetInitialGuess(const real_t*) override;
        virtual void Solve(const real_t, const real_t) override;

        void StoreSolution(len_t iteration);
    };

    PetscErrorCode SNES_convergence_test(SNES, PetscInt, PetscReal, PetscReal, PetscReal, SNESConvergedReason*, void*);
    PetscErrorCode SNES_set_jacobian(SNES, Vec, Mat, Mat, void*);
    PetscErrorCode SNES_set_function(SNES, Vec, Vec, void*);
    void SNES_update_system(const int_t, SolverSNES*);
    PetscErrorCode SNES_solution_obtained(SNES, PetscInt, PetscReal, void*);
}

#endif/*_DREAM_SOLVER_SNES_HPP*/
