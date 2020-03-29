#ifndef _DREAM_FVM_EQUATION_SYSTEM_HPP
#define _DREAM_FVM_EQUATION_SYSTEM_HPP

#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/Equation.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAM {
    class EquationSystem {
    private:
        /// GRIDS
        /// NOTE: These are owned by the parent 'Simulation' object,
        /// and so we should not delete them in the EquationSystem object.
        FVM::RadialGrid *fluidGrid = nullptr;
        FVM::Grid *hottailGrid = nullptr;
        FVM::Grid *runawayGrid = nullptr;

        FVM::BlockMatrix *matrix = nullptr;

    public:
        enum compregion {
            REGION_FLUID,
            REGION_HOTTAIL,
            REGION_RUNAWAY
        };

        EquationSystem(FVM::RadialGrid*, FVM::Grid*, FVM::Grid*);
        ~EquationSystem();

        bool HasHotTailGrid() const { return (this->hottailGrid != nullptr); }
        bool HasRunawayGrid() const { return (this->runawayGrid != nullptr); }

        // Add an unknown to the equation system
        int_t SetUnknown(const std::string&, enum compregion);

        // Set the equation for the specified unknown (blockrow),
        // in the specified block matrix column (blockcol).
        void SetEquation(int blockrow, int blockcol, FVM::Equation&);
    };
}

#endif/*_DREAM_FVM_EQUATION_SYSTEM_HPP*/
