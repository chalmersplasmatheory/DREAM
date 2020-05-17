
#include "FVM/config.h"
#include "CollisionQuantity.hpp"
#include "CoulombLogarithm.hpp"
#include "SlowingDownFrequency.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Constants.hpp"

namespace DREAM {
    class ParallelDiffusionFrequency : public CollisionQuantity{
    private:
        real_t **nonlinearMat = nullptr;
        real_t *trapzWeights = nullptr;
        
        real_t *Tnormalized = nullptr;
        SlowingDownFrequency *nuS;
        CoulombLogarithm *lnLambdaEE;
        real_t rescaleFactor(len_t ir, real_t gamma);
        void calculateIsotropicNonlinearOperatorMatrix();
        void GetNonlinearPartialContribution(const real_t* lnLc, real_t *&partQty);
    protected:
        virtual void AllocatePartialQuantities() override;
        void DeallocatePartialQuantities();        
        virtual void RebuildPlasmaDependentTerms() override;
        virtual void RebuildConstantTerms() override;
        virtual void AssembleQuantity(real_t **&collisionQuantity, len_t nr, len_t np1, len_t np2, len_t fluxGridType) override;

    public:
    
        ParallelDiffusionFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,
            SlowingDownFrequency *nuS, CoulombLogarithm *lnLee,
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantityHandler::collqtyhand_settings *cqset);

        virtual real_t evaluateAtP(len_t ir, real_t p) override;

        void AddNonlinearContribution();


    };
}






