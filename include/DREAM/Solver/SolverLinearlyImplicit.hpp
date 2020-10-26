#ifndef _DREAM_SOLVER_LINEARLY_IMPLICIT_HPP
#define _DREAM_SOLVER_LINEARLY_IMPLICIT_HPP

#include <vector>
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/MatrixInverter.hpp"
#include "FVM/TimeKeeper.hpp"

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
        len_t nTimeStep=0;

        enum OptionConstants::linear_solver linearSolver = OptionConstants::LINEAR_SOLVER_LU;

        FVM::TimeKeeper *timeKeeper;
        len_t timerTot, timerRebuild, timerMatrix, timerInvert;

        // Debug settings
        bool printmatrixinfo = false, savematrix = false, saverhs = false;
        len_t savetimestep = 0;

    protected:
        virtual void initialize_internal(const len_t, std::vector<len_t>&) override;

    public:
        SolverLinearlyImplicit(
            FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*,
            enum OptionConstants::linear_solver ls=OptionConstants::LINEAR_SOLVER_LU
        );
        virtual ~SolverLinearlyImplicit();

        real_t CurrentTime() const { return this->t; }
        real_t CurrentTimeStep() const { return this->dt; }
        FVM::BlockMatrix *GetMatrix() { return this->matrix; }

        virtual void SetInitialGuess(const real_t*) override;
        virtual void Solve(const real_t, const real_t) override;

        virtual void PrintTimings() override;
        virtual void SaveTimings(SFile*, const std::string& path="") override;

        void SaveDebugInfo(len_t, FVM::Matrix*, const real_t*);
        void SetDebugMode(bool, bool, bool, int_t);
    };
}

#endif/*_DREAM_SOLVER_LINEARLY_IMPLICIT_HPP*/
