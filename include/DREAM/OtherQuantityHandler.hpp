#ifndef _OTHER_QUANTITY_HANDLER_HPP
#define _OTHER_QUANTITY_HANDLER_HPP

#include <map>
#include <vector>
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/TransportBC.hpp"
#include "DREAM/OtherQuantity.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/PostProcessor.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/QuantityData.hpp"
#include "DREAM/Settings/OptionConstants.hpp"

#include "DREAM/Equations/Fluid/RadiatedPowerTerm.hpp"
#include "DREAM/Equations/Fluid/OhmicHeatingTerm.hpp"
#include "DREAM/Equations/Fluid/CollisionalEnergyTransferKineticTerm.hpp"
#include "FVM/Equation/AdvectionDiffusionTerm.hpp"

namespace DREAM {
    class OtherQuantityHandler {
    public:
        struct eqn_terms {
            // Terms in the heat equation:
            DREAM::RadiatedPowerTerm *T_cold_radiation=nullptr; 
            DREAM::OhmicHeatingTerm *T_cold_ohmic=nullptr;
            DREAM::CollisionalEnergyTransferKineticTerm *T_cold_fhot_coll=nullptr;
            DREAM::CollisionalEnergyTransferKineticTerm *T_cold_fre_coll=nullptr;
            DREAM::FVM::AdvectionDiffusionTerm *T_cold_transport=nullptr;
            // Radial transport boundary conditions
            DREAM::TransportAdvectiveBC *f_re_advective_bc=nullptr;
            DREAM::TransportDiffusiveBC *f_re_diffusive_bc=nullptr;
            DREAM::TransportAdvectiveBC *f_hot_advective_bc=nullptr;
            DREAM::TransportDiffusiveBC *f_hot_diffusive_bc=nullptr;
            DREAM::TransportAdvectiveBC *n_re_advective_bc=nullptr;
            DREAM::TransportDiffusiveBC *n_re_diffusive_bc=nullptr;
        };

    private:
        std::vector<OtherQuantity*> all_quantities;
        std::vector<OtherQuantity*> registered;

        std::map<std::string, std::vector<std::string>> groups;

        CollisionQuantityHandler *cqtyHottail, *cqtyRunaway;
        PostProcessor *postProcessor;
        RunawayFluid *REFluid;
        FVM::UnknownQuantityHandler *unknowns;
        std::vector<UnknownQuantityEquation*> *unknown_equations;
        FVM::Grid *fluidGrid, *hottailGrid, *runawayGrid, *scalarGrid;

        len_t id_f_hot, id_f_re, id_ncold, id_n_re, id_Tcold, id_Eterm;

        struct eqn_terms *tracked_terms;

    public:
        OtherQuantityHandler(
            CollisionQuantityHandler*, CollisionQuantityHandler*,
            PostProcessor*, RunawayFluid*, FVM::UnknownQuantityHandler*,
            std::vector<UnknownQuantityEquation*>*,
            FVM::Grid*, FVM::Grid*, FVM::Grid*, FVM::Grid*,
            struct eqn_terms*
        );
        ~OtherQuantityHandler();

        void DefineQuantities();
        OtherQuantity *GetByName(const std::string&);
        len_t GetNRegistered() const { return this->registered.size(); }

        bool RegisterGroup(const std::string&);
        void RegisterQuantity(const std::string&);
        void RegisterQuantity(OtherQuantity*);
        void RegisterAllQuantities();
        void StoreAll(const real_t);

        void SaveSFile(SFile*, const std::string& path="other");
    };

    class OtherQuantityException : public DREAM::FVM::FVMException {
    public:
        template<typename ... Args>
        OtherQuantityException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("OtherQuantity");
        }
    };
}

#endif/*_OTHER_QUANTITY_HANDLER_HPP*/
