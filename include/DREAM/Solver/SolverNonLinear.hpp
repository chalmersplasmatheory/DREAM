#ifndef _DREAM_SOLVER_NON_LINEAR_HPP
#define _DREAM_SOLVER_NON_LINEAR_HPP

#include "FVM/config.h"

#include <petsc.h>
#include <vector>
#include "DREAM/ConvergenceChecker.hpp"
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/MatrixInverter.hpp"
#include "FVM/TimeKeeper.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
	class SolverNonLinear : public Solver {
	private:
		FVM::BlockMatrix *jacobian = nullptr;
		FVM::MatrixInverter *inverter = nullptr;
		Vec petsc_F, petsc_dx;

        enum OptionConstants::linear_solver linearSolver = OptionConstants::LINEAR_SOLVER_LU;

		int_t maxiter=100;
		real_t reltol=1e-6;
		bool verbose=false;

		len_t iteration=0, nTimeStep=0;
		real_t t, dt;
		real_t *x0, *x1, *dx;
		real_t *x_2norm, *dx_2norm;

        FVM::TimeKeeper *timeKeeper;
        len_t timerTot, timerRebuild, timerResidual, timerJacobian, timerInvert;

        // Debug settings
        bool printjacobianinfo = false, savejacobian = false, savevector = false, savenumjac = false;
        len_t savetimestep = 0, saveiteration = 1;

	protected:
		virtual void initialize_internal(const len_t, std::vector<len_t>&) override;

        void _EvaluateF(const real_t*, real_t*, FVM::BlockMatrix*);
        void _EvaluateJacobianNumerically(FVM::BlockMatrix*);

	public:
		SolverNonLinear(
			FVM::UnknownQuantityHandler*,
			std::vector<UnknownQuantityEquation*>*,
            enum OptionConstants::linear_solver ls=OptionConstants::LINEAR_SOLVER_LU,
			const int_t maxiter=100, const real_t reltol=1e-6,
			bool verbose=false
		);
		virtual ~SolverNonLinear();

		void Allocate();
		void Deallocate();
		const std::string& GetNonTrivialName(const len_t);

        real_t CurrentTime() const { return this->t; }
        real_t CurrentTimeStep() const { return this->dt; }

		// GETTERS
		len_t GetIteration() const { return this->iteration; }
		int_t MaxIter() const { return this->maxiter; }
		real_t RelTol() const { return this->reltol; }
		bool Verbose() const  { return this->verbose; }

		// Setters
		void SetIteration(const len_t i) { this->iteration = i; }

		bool IsConverged(const real_t*, const real_t*);

		virtual void SetInitialGuess(const real_t*) override;
		virtual void Solve(const real_t, const real_t) override;
		
		void AcceptSolution();
        void SaveNumericalJacobian(const std::string& name="petsc_jacobian");
        // We keep these separate to make it possible to call 'SaveJacobian()'
        // from GDB
        void SaveJacobian();
        void SaveJacobian(const std::string& name);
		void StoreSolution(const real_t*);
		const real_t *TakeNewtonStep();
		const real_t *UpdateSolution(const real_t*);

        virtual void PrintTimings() override;
        virtual void SaveTimings(SFile*, const std::string& path="") override;

        void SaveDebugInfo(len_t, len_t);
        void SetDebugMode(bool, bool, bool, bool, int_t, int_t);
	};
}

#endif/*_DREAM_SOLVER_NON_LINEAR_HPP*/
