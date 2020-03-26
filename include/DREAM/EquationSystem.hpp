#ifndef _DREAM_FVM_EQUATION_SYSTEM_HPP
#define _DREAM_FVM_EQUATION_SYSTEM_HPP

#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/Equation.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAM {
    class EquationSystem {
    private:
        /*********************
         * GRID DEFINITIONS  *
         *********************/
        RadialGrid *fluidGrid;
        Grid *hottailGrid;
        Grid *runawayGrid;

        BlockMatrix *matrix;

    public:
        enum compregion {
            REGION_FLUID,
            REGION_HOTTAIL,
            REGION_RUNAWAY
        };

        EquationSystem();
        ~EquationSystem();

        // Add an unknown to the equation system
        int SetUnknown(const std::string&, enum compregion);

        // Set the equation for the specified unknown (blockrow),
        // in the specified block matrix column (blockcol).
        void SetEquation(int blockrow, int blockcol, Equation&);
    };
}

#endif/*_DREAM_FVM_EQUATION_SYSTEM_HPP*/
