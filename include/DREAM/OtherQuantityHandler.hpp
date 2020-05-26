#ifndef _OTHER_QUANTITY_HANDLER_HPP
#define _OTHER_QUANTITY_HANDLER_HPP

#include <map>
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/QuantityData.hpp"

namespace DREAM {
    class OtherQuantityHandler {
    public:
        enum quantity_id {
            OTHER_QTY_NU_S_HOTTAIL=1,       // Slowing-down collision frequency on hot-tail grid
            OTHER_QTY_NU_S_HOTTAIL_FR=2,    // Slowing-down collision frequency on hot-tail grid
            OTHER_QTY_NU_S_HOTTAIL_F1=3,    // Slowing-down collision frequency on hot-tail grid
            OTHER_QTY_NU_S_HOTTAIL_F2=4,    // Slowing-down collision frequency on hot-tail grid
            OTHER_QTY_NU_S_RUNAWAY=5,       // Slowing-down collision frequency on runaway grid
            OTHER_QTY_NU_S_RUNAWAY_FR=6,    // Slowing-down collision frequency on runaway grid
            OTHER_QTY_NU_S_RUNAWAY_F1=7,    // Slowing-down collision frequency on runaway grid
            OTHER_QTY_NU_S_RUNAWAY_F2=8,    // Slowing-down collision frequency on runaway grid

            // This entry should always come last
            OTHER_QTY_LAST
        };
    private:
        std::map<enum quantity_id, FVM::QuantityData*> data;

        CollisionQuantityHandler *cqtyHottail, *cqtyRunaway;
        FVM::Grid *fluidGrid, *hottailGrid, *runawayGrid;

        FVM::QuantityData *_ConstructQuantity(enum quantity_id);
        const char *_GetName(enum quantity_id);
        void _StoreQuantity(const real_t, enum quantity_id, FVM::QuantityData*);
    public:
        OtherQuantityHandler(
            CollisionQuantityHandler*, CollisionQuantityHandler*,
            FVM::Grid*, FVM::Grid*, FVM::Grid*
        );
        ~OtherQuantityHandler();

        len_t GetNRegistered() const { return this->data.size(); }

        void RegisterQuantity(enum quantity_id);
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
