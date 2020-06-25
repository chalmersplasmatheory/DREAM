#ifndef _DREAM_EQUATION_SCALAR_WALL_CURRENT_TERMS_HPP
#define _DREAM_EQUATION_SCALAR_WALL_CURRENT_TERMS_HPP

#include "FVM/Equation/DiagonalLinearTerm.hpp"
#include "FVM/Equation/ScalarLinearTerm.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {


    /**
     * Implementation of a term which represents the difference
     * between the poloidal flux at r=a and at the wall r=b.
     * It corresponds to 
     * T = I_p(a) * integral(1/(VpVol*FSA_nablaR2OverR2) , r, a, b).
     * TODO: change from the  placeholder value which is only valid in 
     * the cylindrical limit
     */
    class PlasmaEdgeToWallInductanceTerm : public FVM::DiagonalLinearTerm {
    private:
        real_t a, b; // plasma edge radius and wall radius, respectively
    protected:
        virtual void SetWeights() override {
            real_t integralTerm = 1/(4*M_PI*M_PI) * log(b/a);
            weights[0] = -4*M_PI*M_PI*Constants::mu0*integralTerm;
        }
    public:
        PlasmaEdgeToWallInductanceTerm(FVM::Grid* scalarGrid, real_t a, real_t b) 
            : FVM::DiagonalLinearTerm(scalarGrid), a(a), b(b) {}
    };


    /**
     * Implementation of a class which represents the integral
     * I_p = int(...*j_tot,r,0,a)
     * which gives the total toroidal plasma current.
     */
    class TotalPlasmaCurrentFromJTot : public FVM::ScalarLinearTerm {
    protected:
        virtual void SetWeights() override {
            for(len_t i=0; i<nWeights; i++)
                weights[i] = -GetIpIntegrand(i,targetGrid->GetRadialGrid());
        }
    public:
        TotalPlasmaCurrentFromJTot(FVM::Grid* scalarGrid, FVM::Grid* targetGrid,
            FVM::UnknownQuantityHandler *u, const len_t uqtyId) 
            : FVM::ScalarLinearTerm(scalarGrid,targetGrid,u,uqtyId){}
        
        /**
         * Returns the integrand connecting the plasma current I_p to j_tot
         * according to I_p = sum_i j_tot(r_i) integrand_i
         */
        static real_t GetIpIntegrand(const len_t ir, FVM::RadialGrid *rGrid) {
            const real_t dr = rGrid->GetDr(ir);
            const real_t VpVol = rGrid->GetVpVol(ir);
            const real_t Bmin = rGrid->GetBmin(ir);
            const real_t R2inv = rGrid->GetFSA_1OverR2(ir);
            const real_t G = rGrid->GetBTorG(ir);
            return dr*VpVol*G*R2inv/(2*M_PI*Bmin);
        }
    };
}

#endif /*_DREAM_EQUATION_SCALAR_WALL_CURRENT_TERMS_HPP*/
