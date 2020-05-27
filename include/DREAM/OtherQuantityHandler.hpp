#ifndef _OTHER_QUANTITY_HANDLER_HPP
#define _OTHER_QUANTITY_HANDLER_HPP

#include <map>
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/OtherQuantity.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/QuantityData.hpp"

namespace DREAM {
    class OtherQuantityHandler {
    private:
        std::vector<OtherQuantity*> all_quantities;
        std::vector<OtherQuantity*> registered;

        CollisionQuantityHandler *cqtyHottail, *cqtyRunaway;
        FVM::Grid *fluidGrid, *hottailGrid, *runawayGrid;
    public:
        OtherQuantityHandler(
            CollisionQuantityHandler*, CollisionQuantityHandler*,
            FVM::Grid*, FVM::Grid*, FVM::Grid*
        );
        ~OtherQuantityHandler();

        void DefineQuantities();
        OtherQuantity *GetByName(const std::string&);
        len_t GetNRegistered() const { return this->registered.size(); }

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
