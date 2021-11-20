#ifndef _OTHER_QUANTITY_HANDLER_HPP
#define _OTHER_QUANTITY_HANDLER_HPP

namespace DREAM { class OtherQuantityHandler; }

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
#include "DREAM/Equations/Fluid/SvenssonTransport.hpp"
#include "DREAM/Equations/Fluid/CollisionalEnergyTransferREFluidTerm.hpp"
#include "DREAM/Equations/Fluid/HottailRateTerm.hpp"
#include "DREAM/Equations/Fluid/HyperresistiveDiffusionTerm.hpp"
#include "DREAM/Equations/Kinetic/RipplePitchScattering.hpp"
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
            DREAM::CollisionalEnergyTransferREFluidTerm *T_cold_nre_coll=nullptr;
            DREAM::FVM::AdvectionDiffusionTerm *T_cold_transport=nullptr;
            DREAM::FVM::Operator *T_cold_ion_coll=nullptr;
            // Radial transport boundary conditions
            DREAM::TransportAdvectiveBC *f_re_advective_bc=nullptr;
            DREAM::TransportDiffusiveBC *f_re_diffusive_bc=nullptr;
            DREAM::TransportAdvectiveBC *f_hot_advective_bc=nullptr;
            DREAM::TransportDiffusiveBC *f_hot_diffusive_bc=nullptr;
            DREAM::TransportAdvectiveBC *n_re_advective_bc=nullptr;
            DREAM::TransportDiffusiveBC *n_re_diffusive_bc=nullptr;
            DREAM::TransportAdvectiveBC *T_cold_advective_bc=nullptr;
            DREAM::TransportDiffusiveBC *T_cold_diffusive_bc=nullptr;
            // Svensson transport coefficients
            DREAM::SvenssonTransportAdvectionTermA *svensson_A=nullptr;
            DREAM::SvenssonTransportDiffusionTerm  *svensson_D=nullptr;
            DREAM::SvenssonTransportAdvectionTermD *svensson_advD=nullptr;
            // Magnetic ripple pitch scattering
            DREAM::RipplePitchScattering *f_hot_ripple_Dxx=nullptr;
            DREAM::RipplePitchScattering *f_re_ripple_Dxx=nullptr;
            // Runaway rate term
            DREAM::HottailRateTerm *n_re_hottail_rate=nullptr;
            // Hyperresistive diffusion term
            DREAM::HyperresistiveDiffusionTerm *psi_p_hyperresistive=nullptr;
        };

    protected:
        std::vector<OtherQuantity*> all_quantities;
        std::vector<OtherQuantity*> registered;

        std::map<std::string, std::vector<std::string>> groups;

        CollisionQuantityHandler *cqtyHottail, *cqtyRunaway;
        PostProcessor *postProcessor;
        RunawayFluid *REFluid;
        FVM::UnknownQuantityHandler *unknowns;
        std::vector<UnknownQuantityEquation*> *unknown_equations;
        IonHandler *ions;
        FVM::Grid *fluidGrid, *hottailGrid, *runawayGrid, *scalarGrid;

        // indices to unknownquantities
        len_t 
            id_f_hot, id_f_re, id_ncold, id_n_re, id_Tcold, id_Wcold,
            id_Eterm, id_jtot, id_psip, id_Ip, id_psi_edge, id_psi_wall;

        // helper arrays with enough memory allocated to store the hottail and runaway grids 
        real_t *kineticVectorHot; 
        real_t *kineticVectorRE; 

        // helper functions for evaluating other quantities
        real_t integratedKineticBoundaryTerm(
            len_t id_f, std::function<real_t(len_t,len_t,FVM::MomentumGrid*)> momentFunction, FVM::Grid*, 
            FVM::BC::BoundaryCondition*, FVM::BC::BoundaryCondition*, 
            real_t *kineticVector
        );
        real_t evaluateMagneticEnergy();
        real_t integrateWeightedMaxwellian(len_t, real_t, real_t, std::function<real_t(len_t,real_t)>);
        struct eqn_terms *tracked_terms;

    public:
        OtherQuantityHandler(
            CollisionQuantityHandler*, CollisionQuantityHandler*,
            PostProcessor*, RunawayFluid*, FVM::UnknownQuantityHandler*,
            std::vector<UnknownQuantityEquation*>*, IonHandler*,
            FVM::Grid*, FVM::Grid*, FVM::Grid*, FVM::Grid*,
            struct eqn_terms*
        );
        virtual ~OtherQuantityHandler();

        void DefineQuantities();
        OtherQuantity *GetByName(const std::string&);
        len_t GetNRegistered() const { return this->registered.size(); }
        std::vector<OtherQuantity*> GetRegisteredQuantities() { return this->registered; }

        bool RegisterGroup(const std::string&);
        void RegisterQuantity(const std::string&, bool ignorefail=false);
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
