#ifndef _DREAM_EQUATION_SYSTEM_HPP
#define _DREAM_EQUATION_SYSTEM_HPP

namespace DREAM { class EquationSystem; }

#include <map>
#include <string>
#include <vector>
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/OtherQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/TimeStepper/TimeStepper.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/Equation.hpp"
#include "FVM/FVMException.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/UnknownQuantity.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "IonHandler.hpp"
#include "FVM/QuantityData.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"

namespace DREAM {
    class EquationSystem {
    private:
        /// GRIDS
        /// NOTE: These are owned by the parent 'Simulation' object,
        /// and so we should not delete them in the EquationSystem object.
        FVM::Grid *fluidGrid = nullptr;
        FVM::Grid *hottailGrid = nullptr;
        FVM::Grid *runawayGrid = nullptr;

        enum OptionConstants::momentumgrid_type hottailGrid_type;
        enum OptionConstants::momentumgrid_type runawayGrid_type;

        IonHandler *ionHandler=nullptr;
        Solver *solver=nullptr;
        TimeStepper *timestepper=nullptr;

        FVM::UnknownQuantityHandler unknowns;
        std::vector<UnknownQuantityEquation*> unknown_equations;
        std::vector<len_t> nontrivial_unknowns;

        CollisionQuantityHandler *cqh_hottail=nullptr;
        CollisionQuantityHandler *cqh_runaway=nullptr;

        RunawayFluid *REFluid = nullptr;

        OtherQuantityHandler *otherQuantityHandler=nullptr;

        real_t currentTime;
        std::vector<real_t> times;

        len_t matrix_size=0;

    public:
        EquationSystem(FVM::Grid*, enum OptionConstants::momentumgrid_type, FVM::Grid*, enum OptionConstants::momentumgrid_type, FVM::Grid*);
        ~EquationSystem();

        FVM::Grid *GetFluidGrid() { return this->fluidGrid; }
        FVM::Grid *GetHotTailGrid() { return this->hottailGrid; }
        FVM::Grid *GetRunawayGrid() { return this->runawayGrid; }

        enum OptionConstants::momentumgrid_type GetHotTailGridType()
        { return this->hottailGrid_type; }
        enum OptionConstants::momentumgrid_type GetRunawayGridType()
        { return this->runawayGrid_type; }

        bool HasInitialValue(const len_t unknId) const { return this->unknowns.HasInitialValue(unknId); }
        bool HasHotTailGrid() const { return (this->hottailGrid != nullptr); }
        bool HasRunawayGrid() const { return (this->runawayGrid != nullptr); }

        CollisionQuantityHandler *GetHotTailCollisionHandler() { return this->cqh_hottail; }
        CollisionQuantityHandler *GetRunawayCollisionHandler() { return this->cqh_runaway; }
        
        RunawayFluid *GetREFluid() { return this->REFluid; }

        OtherQuantityHandler *GetOtherQuantityHandler() { return this->otherQuantityHandler; }

        FVM::UnknownQuantity *GetUnknown(const len_t i) { return unknowns.GetUnknown(i); }
        FVM::UnknownQuantityHandler *GetUnknownHandler() { return &unknowns; }
        IonHandler *GetIonHandler() { return this->ionHandler; }
        std::vector<len_t> *GetNonTrivialUnknowns() { return &nontrivial_unknowns; }
        UnknownQuantityEquation *GetEquation(const len_t i) { return unknown_equations.at(i); }
        std::vector<UnknownQuantityEquation*> *GetEquations() { return &unknown_equations; }

        real_t *GetUnknownData(const len_t i) { return unknowns.GetUnknownData(i); }
        len_t GetUnknownID(const std::string& name) { return unknowns.GetUnknownID(name); }
        len_t GetNUnknowns() const { return this->unknowns.GetNUnknowns(); }

        void ProcessSystem();

        // Add an unknown to the equation system
        len_t SetUnknown(const std::string& name, FVM::Grid *grid, const len_t nMultiples=1)
        { return unknowns.InsertUnknown(name, grid, nMultiples); }

        // Set the equation for the specified unknown (blockrow),
        // in the specified block matrix column (blockcol).
        void SetEquation(len_t blockrow, len_t blockcol, FVM::Equation *eqn, const std::string& desc="");
        //{ return unknowns.SetEquation(blockrow, blockcol, eqn); }

        // Set equation by name of the unknown
        // NOTE: These are slower and should be used only when
        // performance is not a concern
        void SetEquation(len_t blockrow, const std::string&, FVM::Equation*, const std::string& desc="");
        void SetEquation(const std::string&, len_t blockcol, FVM::Equation*, const std::string& desc="");
        void SetEquation(const std::string&, const std::string&, FVM::Equation*, const std::string& desc="");

        void SetHotTailCollisionHandler(CollisionQuantityHandler *cqh)
        { this->cqh_hottail = cqh; }
        void SetRunawayCollisionHandler(CollisionQuantityHandler *cqh)
        { this->cqh_runaway = cqh; }

        void SetREFluid(RunawayFluid *REF)
        { this->REFluid = REF; }

        void SetInitialValue(const len_t, const real_t*, const real_t t0=0);
        void SetInitialValue(const std::string&, const real_t*, const real_t t0=0);

        void SetIonHandler(IonHandler *ih) { this->ionHandler = ih; }
        void SetOtherQuantityHandler(OtherQuantityHandler *oqh) { this->otherQuantityHandler = oqh; }
        void SetSolver(Solver*);
        void SetTimeStepper(TimeStepper *ts) { this->timestepper = ts; }

        void SaveGrids(SFile*, const std::string& path="");
        void SaveIonMetaData(SFile*, const std::string& path="");
        void SaveMomentumGrid(SFile*, const std::string&, FVM::Grid*, enum OptionConstants::momentumgrid_type);
        void WriteCopyArray(SFile*, const std::string&, const real_t *const*, const len_t, const len_t);

        void Solve();
        void Rebuild();
        // Info routines
        void PrintNonTrivialUnknowns();
        void PrintTrivialUnknowns();
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

#endif/*_DREAM_EQUATION_SYSTEM_HPP*/
