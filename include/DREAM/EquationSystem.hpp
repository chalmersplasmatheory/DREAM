#ifndef _DREAM_FVM_EQUATION_SYSTEM_HPP
#define _DREAM_FVM_EQUATION_SYSTEM_HPP

#include <map>
#include <string>
#include <vector>
#include "DREAM/QuantityData.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/Equation.hpp"
#include "FVM/FVMException.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAM {
    class EquationSystem {
    public:
        struct unknown_qty {
            // Name of quantity
            std::string name;
            // Pointer to grid on which the quantity is defined
            FVM::Grid *grid;
            // List of equations associated with the quantity
            std::map<len_t, FVM::Equation*> equations;
            // (Solution) data handler
            QuantityData *data;

            unknown_qty(const std::string& name, FVM::Grid *grid) {
                this->name = name;
                this->grid = grid;
                this->data = new QuantityData(grid);
            }
            ~unknown_qty() {
                delete data;
            }
        };
    private:
        /// GRIDS
        /// NOTE: These are owned by the parent 'Simulation' object,
        /// and so we should not delete them in the EquationSystem object.
        FVM::Grid *fluidGrid = nullptr;
        FVM::Grid *hottailGrid = nullptr;
        FVM::Grid *runawayGrid = nullptr;

        FVM::BlockMatrix *matrix = nullptr;

        std::vector<unknown_qty*> unknowns;

    public:
        EquationSystem(FVM::Grid*, FVM::Grid*, FVM::Grid*);
        ~EquationSystem();

        FVM::Grid *GetFluidGrid() { return this->fluidGrid; }
        FVM::Grid *GetHotTailGrid() { return this->hottailGrid; }
        FVM::Grid *GetRunawayGrid() { return this->runawayGrid; }

        bool HasHotTailGrid() const { return (this->hottailGrid != nullptr); }
        bool HasRunawayGrid() const { return (this->runawayGrid != nullptr); }

        real_t *GetUnknownData(const len_t);
        len_t GetUnknownID(const std::string&);

        // Add an unknown to the equation system
        len_t SetUnknown(const std::string&, FVM::Grid*);

        // Set the equation for the specified unknown (blockrow),
        // in the specified block matrix column (blockcol).
        void SetEquation(len_t blockrow, len_t blockcol, FVM::Equation&);
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
