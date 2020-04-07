#ifndef _DREAM_FVM_EQUATION_SYSTEM_HPP
#define _DREAM_FVM_EQUATION_SYSTEM_HPP

#include <map>
#include <string>
#include <vector>
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "DREAM/TimeStepper/TimeStepper.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/Equation.hpp"
#include "FVM/FVMException.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/UnknownQuantity.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/QuantityData.hpp"

namespace DREAM {
    class EquationSystem {
    private:
        /// GRIDS
        /// NOTE: These are owned by the parent 'Simulation' object,
        /// and so we should not delete them in the EquationSystem object.
        FVM::Grid *fluidGrid = nullptr;
        FVM::Grid *hottailGrid = nullptr;
        FVM::Grid *runawayGrid = nullptr;

        Solver *solver;
        TimeStepper *timestepper;

        FVM::UnknownQuantityHandler unknowns;
        std::vector<UnknownQuantityEquation*> unknown_equations;
        std::vector<len_t> nontrivial_unknowns;

        real_t currentTime;

    public:
        EquationSystem(FVM::Grid*, FVM::Grid*, FVM::Grid*);
        ~EquationSystem();

        FVM::Grid *GetFluidGrid() { return this->fluidGrid; }
        FVM::Grid *GetHotTailGrid() { return this->hottailGrid; }
        FVM::Grid *GetRunawayGrid() { return this->runawayGrid; }

        bool HasHotTailGrid() const { return (this->hottailGrid != nullptr); }
        bool HasRunawayGrid() const { return (this->runawayGrid != nullptr); }

        FVM::UnknownQuantity *GetUnknown(const len_t i) { return unknowns.GetUnknown(i); }
        FVM::UnknownQuantityHandler *GetUnknownHandler() { return &unknowns; }
        std::vector<len_t> *GetNonTrivialUnknowns() { return &nontrivial_unknowns; }
        std::vector<UnknownQuantityEquation*> *GetEquations() { return &unknown_equations; }

        real_t *GetUnknownData(const len_t i) { return unknowns.GetUnknownData(i); }
        len_t GetUnknownID(const std::string& name) { return unknowns.GetUnknownID(name); }
        len_t GetNUnknowns() const { return this->unknowns.GetNUnknowns(); }

        void ProcessSystem();

        // Add an unknown to the equation system
        len_t SetUnknown(const std::string& name, FVM::Grid *grid)
        { return unknowns.InsertUnknown(name, grid); }

        // Set the equation for the specified unknown (blockrow),
        // in the specified block matrix column (blockcol).
        void SetEquation(len_t blockrow, len_t blockcol, FVM::Equation *eqn);
        //{ return unknowns.SetEquation(blockrow, blockcol, eqn); }

        // Set equation by name of the unknown
        // NOTE: These are slower and should be used only when
        // performance is not a concern
        void SetEquation(len_t blockrow, const std::string&, FVM::Equation*);
        void SetEquation(const std::string&, len_t blockcol, FVM::Equation*);
        void SetEquation(const std::string&, const std::string&, FVM::Equation*);

        void SetSolver(Solver *solver) { this->solver = solver; }
        void SetTimeStepper(TimeStepper *ts) { this->timestepper = ts; }

        void Solve();
    };

    class EquationSystemException : public DREAM::FVM::FVMException {
    public:
        template<typename ... Args>
        EquationSystemException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("EquationSystem");
        }
    };
}

#endif/*_DREAM_FVM_EQUATION_SYSTEM_HPP*/
