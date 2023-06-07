#ifndef _DREAM_SOLVER_LINEARLY_IMPLICIT_HPP
#define _DREAM_SOLVER_LINEARLY_IMPLICIT_HPP

#include <vector>
#include "DREAM/EquationSystem.hpp"
#include "DREAM/OutputGenerator.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/BlockMatrix.hpp"
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

        EquationSystem *eqsys;

        real_t t, dt;
        len_t nTimeStep=0;

        FVM::TimeKeeper *timeKeeper;
        len_t timerTot, timerRebuild, timerMatrix, timerInvert;

		// Iterations of external iterator
		std::vector<len_t> extiter_nIterations;
		
        // Debug settings
        bool printmatrixinfo = false, savematrix = false, saverhs = false, savesystem = false;
        len_t savetimestep = 0;

    protected:
        virtual void initialize_internal(const len_t, std::vector<len_t>&) override;

    public:
        SolverLinearlyImplicit(
            FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*,
            EquationSystem*, const bool verbose=false,
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
        void SetDebugMode(bool, bool, bool, int_t, bool);

        virtual void WriteDataSFile(SFile*, const std::string&) override;
    };
}

#endif/*_DREAM_SOLVER_LINEARLY_IMPLICIT_HPP*/
