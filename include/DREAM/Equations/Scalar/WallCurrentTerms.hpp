#ifndef _DREAM_EQUATION_SCALAR_WALL_CURRENT_TERMS_HPP
#define _DREAM_EQUATION_SCALAR_WALL_CURRENT_TERMS_HPP

#include "FVM/Equation/ScalarLinearTerm.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

/**
 * Implementation of a class which represents the j_||/(B/Bmin) term in Ampere's law.
 */
namespace DREAM {
    class PoloidalFluxAtEdgeTerm : public FVM::ScalarLinearTerm {
    protected:
        virtual void SetWeights() override {weights[nWeights-1] = -1;}
    public:
        PoloidalFluxAtEdgeTerm(FVM::Grid* scalarGrid, FVM::Grid* targetGrid,
            FVM::UnknownQuantityHandler *u, const len_t uqtyId) 
            : FVM::ScalarLinearTerm(scalarGrid,targetGrid,u,uqtyId){}
    };
}

namespace DREAM {
    class SOLMutualInductanceTerm : public FVM::ScalarLinearTerm {
    private:
        real_t a,b; // plasma and wall radius, respectively
    protected:
        virtual void SetWeights() override {
            real_t integralTerm = 1/(4*M_PI*M_PI) * log(b/a);
            // integralterm should be the radial integral of
            // 1/(VpVol*FSA_nablaR2OverR2) from a to b.
            // Now set to the cylindrical limit as a placeholder.
            weights[0] = -2*M_PI*Constants::mu0*integralTerm;
        }
    public:
        SOLMutualInductanceTerm(FVM::Grid* scalarGrid, FVM::Grid* targetGrid,
            FVM::UnknownQuantityHandler *u, const len_t uqtyId,
            real_t a, real_t b) 
            : FVM::ScalarLinearTerm(scalarGrid,targetGrid,u,uqtyId),
              a(a), b(b) {}

    };
}


namespace DREAM {
    class TotalPlasmaCurrentFromJTot : public FVM::ScalarLinearTerm {
    protected:
        virtual void SetWeights() override {
            FVM::RadialGrid *rGrid = targetGrid->GetRadialGrid();
            const real_t *dr = rGrid->GetDr();
            for(len_t i=0; i<nWeights; i++){
                weights[i] = -dr[i]*rGrid->GetVpVol(i)
                    *rGrid->GetBTorG(i)/rGrid->GetBmin(i)
                    *rGrid->GetFSA_1OverR2(i);
            }
        }
    public:
        TotalPlasmaCurrentFromJTot(FVM::Grid* scalarGrid, FVM::Grid* targetGrid,
            FVM::UnknownQuantityHandler *u, const len_t uqtyId) 
            : FVM::ScalarLinearTerm(scalarGrid,targetGrid,u,uqtyId){}
    };
}

#endif /*_DREAM_EQUATION_SCALAR_WALL_CURRENT_TERMS_HPP*/
