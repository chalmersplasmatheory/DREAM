
#include "FVM/config.h"
#include "CollisionQuantity.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Constants.hpp"
#include "SlowingDownFrequency.hpp"

namespace DREAM {
    class ParallelDiffusionFrequency : public CollisionQuantity{
    private:
        SlowingDownFrequency *nuS;
        void rescaleFrequency(len_t id_unknown, len_t ir, real_t p, real_t *&partQty);
        real_t rescaleFactor(len_t ir, real_t p);
        void calculateIsotropicNonlinearOperatorMatrix();
    protected:
        virtual void GetPartialContribution(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) override;
        virtual void GetPartialContribution_fr(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) override;
        virtual void GetPartialContribution_f1(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) override;
        virtual void GetPartialContribution_f2(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) override;

        virtual void RebuildConstantTerms() override;


        virtual void AllocatePartialQuantities() override {return;}
        virtual void RebuildPlasmaDependentTerms() override {return;}
    public:
    
        ParallelDiffusionFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,
            SlowingDownFrequency *nuS,  
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantityHandler::collqtyhand_settings *cqset);

        virtual real_t evaluateAtP(len_t ir, real_t p) override;

    };
}






