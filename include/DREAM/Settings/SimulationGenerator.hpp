#ifndef _DREAM_PROCESS_SETTINGS_HPP
#define _DREAM_PROCESS_SETTINGS_HPP

#include "DREAM/ADAS.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/IonInterpolator1D.hpp"
#include "DREAM/NIST.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Simulation.hpp"
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/Solver/SolverLinearlyImplicit.hpp"
#include "DREAM/Solver/SolverSNES.hpp"
#include "DREAM/TimeStepper/TimeStepper.hpp"
#include "DREAM/TimeStepper/TimeStepperConstant.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/Interpolator3D.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"

namespace DREAM {
    class SimulationGenerator {
    public:
        // PUBLIC INTERFACE
        static Settings *CreateSettings() {
            Settings *s = new Settings();
            DefineOptions(s);

            return s;
        }
        static void DefineOptions(Settings*);
        static Simulation *ProcessSettings(Settings*);

        // FOR INTERNAL USE
        static EquationSystem *ConstructEquationSystem(Settings*, FVM::Grid*, FVM::Grid*,  enum OptionConstants::momentumgrid_type, FVM::Grid*, enum OptionConstants::momentumgrid_type, FVM::Grid*, ADAS*, NIST*);
        static FVM::Grid *ConstructScalarGrid();
        static FVM::Grid *ConstructHotTailGrid(Settings*, FVM::RadialGrid*, enum OptionConstants::momentumgrid_type*);
        static FVM::Grid *ConstructRunawayGrid(Settings*, FVM::RadialGrid*, FVM::Grid*, enum OptionConstants::momentumgrid_type*);
        
        static RunawayFluid *ConstructRunawayFluid(FVM::Grid *g,
                FVM::UnknownQuantityHandler *unknowns, IonHandler *ih, 
                OptionConstants::momentumgrid_type gridtype, Settings *s);

        static FVM::Grid *ConstructRadialGrid(Settings*);
        static FVM::RadialGrid *ConstructRadialGrid_Cylindrical(const int_t, Settings*);

        static FVM::PXiGrid::PXiMomentumGrid *Construct_PXiGrid(
            Settings*, const std::string&, const real_t, FVM::RadialGrid*
        );

        static enum FVM::Interpolator3D::momentumgrid_type GetInterp3DMomentumGridType(
            enum OptionConstants::momentumgrid_type
        );

        static void DefineOptions_RadialGrid(Settings*);

        static void DefineOptions_ADAS(Settings*);
        static void DefineOptions_CollisionQuantityHandler(Settings*);
        
        static void DefineOptions_RunawayFluid(Settings*);
        static void DefineOptions_EquationSystem(Settings*);
        static void DefineOptions_f_hot(Settings*);
        static void DefineOptions_ElectricField(Settings*);
        static void DefineOptions_T_cold(Settings*);
        static void DefineOptions_j_tot(Settings*);
        //static void DefineOptions_psi_p(Settings*);
        static void DefineOptions_HotTailGrid(Settings*);
        static void DefineOptions_Ions(Settings*);
        static void DefineOptions_KineticGrid(const std::string&, Settings*);
        static void DefineOptions_OtherQuantities(Settings*);
        static void DefineOptions_RunawayGrid(Settings*);
        static void DefineOptions_Solver(Settings*);
        static void DefineOptions_TimeStepper(Settings*);

        static ADAS *LoadADAS(Settings*);
        static NIST *LoadNIST(Settings*);
        static CollisionQuantityHandler *ConstructCollisionQuantityHandler(enum OptionConstants::momentumgrid_type, FVM::Grid *,FVM::UnknownQuantityHandler *, IonHandler *,  Settings*);
        static void ConstructEquations(EquationSystem*, Settings*, ADAS*, NIST*);
        static void ConstructOtherQuantityHandler(EquationSystem*, Settings*);
        static void ConstructSolver(EquationSystem*, Settings*);
        static void ConstructTimeStepper(EquationSystem*, Settings*);
        static void ConstructUnknowns(EquationSystem*, Settings*, FVM::Grid*, FVM::Grid*, FVM::Grid*);
        
        // Routines for constructing specific equations
        static void ConstructEquation_E_field(EquationSystem*, Settings*);
        static void ConstructEquation_E_field_prescribed(EquationSystem*, Settings*);
        static void ConstructEquation_E_field_selfconsistent(EquationSystem*, Settings*);

        static void ConstructEquation_f_hot(EquationSystem*, Settings*);
        static void ConstructEquation_f_hot_maxwellian(EquationSystem*, FVM::Grid*, const real_t*, const real_t*);

        static void ConstructEquation_Ions(EquationSystem*, Settings*, ADAS*);

        static void ConstructEquation_n_cold(EquationSystem*, Settings*);
        static void ConstructEquation_n_cold_prescribed(EquationSystem*, Settings*);
        static void ConstructEquation_n_cold_selfconsistent(EquationSystem*, Settings*);

        static void ConstructEquation_n_hot(EquationSystem*, Settings*);
        static void ConstructEquation_j_hot(EquationSystem*, Settings*);
        static void ConstructEquation_j_ohm(EquationSystem*, Settings*);
        static void ConstructEquation_j_tot(EquationSystem*, Settings*);

        static void ConstructEquation_psi_p(EquationSystem*, Settings*);
        static void ConstructEquations_I_wall(EquationSystem*, Settings*);



        static void ConstructEquation_n_re(EquationSystem*, Settings*);

        static void ConstructEquation_n_tot(EquationSystem*, Settings*);

        static void ConstructEquation_T_cold(EquationSystem*, Settings*, ADAS*, NIST*);
        static void ConstructEquation_T_cold_prescribed(EquationSystem*, Settings*);
        static void ConstructEquation_T_cold_selfconsistent(EquationSystem*, Settings*, ADAS*, NIST*);
        static void ConstructEquation_W_cold(EquationSystem*, Settings*, NIST*);

        // Routines for constructing time steppers
        static TimeStepperConstant *ConstructTimeStepper_constant(Settings*, FVM::UnknownQuantityHandler*);

        // Data loading routines
        static void DefineDataR(const std::string&, Settings*, const std::string& name="data");
        static real_t *LoadDataR(const std::string&, FVM::RadialGrid*, Settings*, const std::string& name="data");
        static void DefineDataRT(const std::string&, Settings*, const std::string& name="data");
        static FVM::Interpolator1D *LoadDataRT(const std::string&, FVM::RadialGrid*, Settings*, const std::string& name="data");
        static void DefineDataR2P(const std::string&, Settings*, const std::string& name="data");
        static FVM::Interpolator3D *LoadDataR2P(const std::string&, Settings*, const std::string& name="data");
        static void DefineDataIonR(const std::string&, Settings*, const std::string& name="data");
        static real_t *LoadDataIonR(const std::string&, FVM::RadialGrid*, Settings*, const len_t, const std::string& name="data");
        static void DefineDataIonRT(const std::string&, Settings*, const std::string& name="data");
        static IonInterpolator1D *LoadDataIonRT(const std::string&, FVM::RadialGrid*, Settings*, const len_t, const std::string& name="data");

        static len_t GetNumberOfIonChargeStates(Settings*);

        // Routines for constructing solvers
        static SolverLinearlyImplicit *ConstructSolver_linearly_implicit(Settings*, FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*);
        static SolverSNES *ConstructSolver_nonlinear_snes(Settings*, FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*);
    };
}

#endif/*_DREAM_PROCESS_SETTINGS_HPP*/
