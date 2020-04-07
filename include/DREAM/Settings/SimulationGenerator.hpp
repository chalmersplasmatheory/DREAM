#ifndef _DREAM_PROCESS_SETTINGS_HPP
#define _DREAM_PROCESS_SETTINGS_HPP

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Simulation.hpp"
#include "DREAM/TimeStepper/TimeStepper.hpp"
#include "DREAM/TimeStepper/TimeStepperConstant.hpp"
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/Solver/SolverLinearlyImplicit.hpp"
#include "DREAM/Solver/SolverSNES.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAM {
    class SimulationGenerator {
    public:
        #include "SimulationGenerator.enum.hpp"

        // PUBLIC INTERFACE
        static Settings *CreateSettings() {
            Settings *s = new Settings();
            DefineOptions(s);

            return s;
        }
        static void DefineOptions(Settings*);
        static Simulation *ProcessSettings(Settings*);

        // FOR INTERNAL USE
        static EquationSystem *ConstructEquationSystem(Settings*, FVM::Grid*, FVM::Grid*, FVM::Grid*);
        static FVM::Grid *ConstructHotTailGrid(Settings*, FVM::RadialGrid*);
        static FVM::Grid *ConstructRunawayGrid(Settings*, FVM::RadialGrid*, FVM::Grid*);

        static FVM::Grid *ConstructRadialGrid(Settings*);
        static FVM::RadialGrid *ConstructRadialGrid_Cylindrical(const int_t, Settings*);

        static FVM::PXiGrid::PXiMomentumGrid *Construct_PXiGrid(
            Settings*, const std::string&, const real_t, FVM::RadialGrid*
        );

        static void DefineOptions_RadialGrid(Settings*);

        static void DefineOptions_EquationSystem(Settings*);
        static void DefineOptions_KineticGrid(const std::string&, Settings*);
        static void DefineOptions_HotTailGrid(Settings*);
        static void DefineOptions_RunawayGrid(Settings*);
        static void DefineOptions_Solver(Settings*);
        static void DefineOptions_TimeStepper(Settings*);

        static void ConstructEquations(EquationSystem*, Settings*);
        static void ConstructSolver(EquationSystem*, Settings*);
        static void ConstructTimeStepper(EquationSystem*, Settings*);
        static void ConstructUnknowns(EquationSystem*, Settings*, FVM::Grid*, FVM::Grid*, FVM::Grid*);

        // Routines for constructing specific equations
        static void ConstructEquation_n_cold(EquationSystem*, Settings*);
        static void ConstructEquation_n_cold_prescribed(EquationSystem*, Settings*);

        // Routines for constructing time steppers
        static TimeStepperConstant *ConstructTimeStepper_constant(Settings*, FVM::UnknownQuantityHandler*);

        // Routines for constructing solvers
        static SolverLinearlyImplicit *ConstructSolver_linearly_implicit(Settings*, FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*);
        static SolverSNES *ConstructSolver_nonlinear_snes(Settings*, FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*);

        // CONSTANTS
        static const char *UQTY_E_FIELD;
        static const char *UQTY_F_HOT;
        static const char *UQTY_F_RE;
        static const char *UQTY_T_COLD;
        static const char *UQTY_N_COLD;
        static const char *UQTY_N_HOT;
        static const char *UQTY_N_RE;
        static const char *UQTY_ION_SPECIES;
        
    };
}

#endif/*_DREAM_PROCESS_SETTINGS_HPP*/
