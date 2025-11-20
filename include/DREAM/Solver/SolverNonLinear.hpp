#ifndef _DREAM_SOLVER_NON_LINEAR_HPP
#define _DREAM_SOLVER_NON_LINEAR_HPP

#include "FVM/config.h"

#include "DREAM/ConvergenceChecker.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Solver/NewtonStepAdjuster.hpp"
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/TimeKeeper.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include <petsc.h>
#include <string>
#include <vector>

namespace DREAM {
	class SolverNonLinear : public Solver {
	private:
		FVM::BlockMatrix *jacobian = nullptr;
		Vec petsc_F, petsc_dx;

		len_t maxiter=100;
		real_t reltol=1e-6;

		len_t iteration=0, nTimeStep=0;
		real_t t, dt;
		real_t *x0, *x1, *dx, *xinit;
		real_t *x_2norm, *dx_2norm;

        FVM::TimeKeeper *timeKeeper;
        len_t timerTot, timerRebuild, timerResidual, timerJacobian, timerInvert;

		bool checkResidual = true;

        NewtonStepAdjuster *adjuster;
        enum OptionConstants::newton_step_adjuster stepAdjusterType = OptionConstants::NEWTON_STEP_ADJUSTER_PHYSICAL;
		std::vector<std::string> stepAdjusterMonitor;

        // Debug settings
        bool printjacobianinfo = false, savejacobian = false, savesolution = false,
            savevector = false, savenumjac = false, savesystem = false, debugrescaled = false;
        len_t savetimestep = 0, saveiteration = 1;

        std::vector<real_t> solver_time;
        std::vector<len_t> nIterations;
        std::vector<bool> usedBackupInverter;

	protected:
		virtual void initialize_internal(const len_t, std::vector<len_t>&) override;

        void _EvaluateF(const real_t*, real_t*, FVM::BlockMatrix*);
        void _EvaluateJacobianNumerically(FVM::BlockMatrix*);
        void _InternalSolve();

	public:
		SolverNonLinear(
			FVM::UnknownQuantityHandler*,
			std::vector<UnknownQuantityEquation*>*, EquationSystem*,
            enum OptionConstants::linear_solver ls=OptionConstants::LINEAR_SOLVER_LU,
            enum OptionConstants::linear_solver bk=OptionConstants::LINEAR_SOLVER_NONE,
            enum OptionConstants::newton_step_adjuster nsa=OptionConstants::NEWTON_STEP_ADJUSTER_PHYSICAL,
			const std::vector<std::string> &nsaMonitor=std::vector<std::string>(),
			const int_t maxiter=100, const real_t reltol=1e-6,
			bool verbose=false, bool checkResidual=true
		);
		virtual ~SolverNonLinear();

		void Allocate();
        void AllocateJacobianMatrix();
		void Deallocate();
		const std::string& GetNonTrivialName(const len_t);
        void InitStepAdjuster(enum OptionConstants::newton_step_adjuster, std::vector<std::string>&);

        real_t CurrentTime() const { return this->t; }
        real_t CurrentTimeStep() const { return this->dt; }

		// GETTERS
		len_t GetIteration() const { return this->iteration; }
		len_t MaxIter() const { return this->maxiter; }
		real_t RelTol() const { return this->reltol; }

		// Setters
		void SetIteration(const len_t i) { this->iteration = i; }

		bool IsConverged(const real_t*, const real_t*);
		bool IsResidualConverged();

		virtual void SetInitialGuess(const real_t*) override;
		virtual void Solve(const real_t, const real_t) override;
        void ResetSolution();
		
		void AcceptSolution();
        void SaveNumericalJacobian(const std::string& name="petsc_jacobian");
        // We keep these separate to make it possible to call 'SaveJacobian()'
        // from GDB
        void SaveJacobian();
        void SaveJacobian(const std::string& name);
		void StoreSolution(const real_t*);
		void EvaluateResidual();
		const real_t *EvaluateNewtonStep();

        virtual void PrintTimings() override;
        virtual void SaveTimings(SFile*, const std::string& path="") override;

        void SaveDebugInfoBefore(len_t, len_t);
        void SaveDebugInfoAfter(len_t, len_t);
        void SetDebugMode(bool, bool, bool, bool, bool, int_t, int_t, bool, bool);

        // Routines for backtracking
        void EvaluateTargetFunction(Vec*);
        void PushTargetFunction(const real_t);

        virtual void SwitchToBackupInverter() override;

        virtual void WriteDataSFile(SFile*, const std::string&) override;
	};
}

#endif/*_DREAM_SOLVER_NON_LINEAR_HPP*/
