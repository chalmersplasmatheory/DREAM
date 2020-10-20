#ifndef _OTHER_QUANTITY_HANDLER_HPP
#define _OTHER_QUANTITY_HANDLER_HPP

#include <map>
#include <vector>
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
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

namespace DREAM {
    class OtherQuantityHandler {
    public:
        struct eqn_terms {
            // Radiated power term in self-consistent T_cold
            // Terms in the heat equation:
            DREAM::RadiatedPowerTerm *T_cold_radiation=nullptr; 
            DREAM::OhmicHeatingTerm *T_cold_ohmic=nullptr;
            DREAM::CollisionalEnergyTransferKineticTerm *T_cold_fhot_coll=nullptr;
            DREAM::CollisionalEnergyTransferKineticTerm *T_cold_fre_coll=nullptr;
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
        FVM::Grid *fluidGrid, *hottailGrid, *runawayGrid;

        len_t id_Eterm;
        len_t id_ncold;
        len_t id_Tcold;

        struct eqn_terms *tracked_terms;

    public:
        OtherQuantityHandler(
            CollisionQuantityHandler*, CollisionQuantityHandler*,
            PostProcessor*, RunawayFluid*, FVM::UnknownQuantityHandler*,
            std::vector<UnknownQuantityEquation*>*,
            FVM::Grid*, FVM::Grid*, FVM::Grid*,
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
