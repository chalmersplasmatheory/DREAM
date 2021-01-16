#ifndef _DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_HPP
#define _DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_HPP

namespace DREAM { class AnalyticDistribution; }

#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Grid/fluxGridType.enum.hpp"
#include "DREAM/Constants.hpp"

/**
 * Base class for analytic distribution functions,
 * containing methods to evaluate distributions
 * at radial index ir, momentum p and pitch xi0.
 * If the df/dx arguments are provided, the methods
 * also return the derivative of the distribution
 * with respect to x.  
 */

namespace DREAM {
    class AnalyticDistribution {
        protected:
            FVM::RadialGrid *rGrid;
            FVM::UnknownQuantityHandler *unknowns;

        public:
            AnalyticDistribution(FVM::RadialGrid *rGrid, FVM::UnknownQuantityHandler *u) 
                : rGrid(rGrid), unknowns(u) {}

            virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*){}
            virtual bool GridRebuilt(){return false;}
            /**
             * Returns the full distribution function F such that
             *  int(Vp*F dp dxi0 dr) 
             * is the total particle number
             */
            virtual real_t evaluateFullDistribution(len_t ir, real_t xi0, real_t p, real_t *dfdxi0=nullptr, real_t *dfdp=nullptr, real_t *dfdr=nullptr) = 0;
            virtual real_t evaluateFullDistribution_fr(len_t ir, real_t xi0, real_t p, real_t *dfdxi0=nullptr, real_t *dfdp=nullptr, real_t *dfdr=nullptr) = 0;
            virtual real_t evaluateFullDistribution(len_t ir, real_t xi0, real_t p, FVM::fluxGridType fgType, real_t *dfdxi0=nullptr, real_t *dfdp=nullptr, real_t *dfdr=nullptr) {
                if (fgType==FVM::FLUXGRIDTYPE_RADIAL)
                    return evaluateFullDistribution_fr(ir,xi0,p,dfdxi0, dfdp, dfdr);
                else 
                    return evaluateFullDistribution(ir,xi0,p,dfdxi0, dfdp, dfdr);
            }

            /** 
             * Returns the pitch-angle averaged distribution function F0
             * evaluated at radial index ir and momentum p, normalized 
             * such that int(4*pi*p^2*F0 dp,0,inf) is the flux-surface 
             * averaged density.
             */
            virtual real_t evaluateEnergyDistribution(len_t ir, real_t p, real_t *dfdp=nullptr, real_t *dfdr=nullptr) = 0;
            virtual real_t evaluateEnergyDistribution_fr(len_t ir, real_t p, real_t *dfdp=nullptr, real_t *dfdr=nullptr) = 0;
            virtual real_t evaluateEnergyDistribution(len_t ir, real_t p, FVM::fluxGridType fgType, real_t *dfdp=nullptr, real_t *dfdr=nullptr) {
                if (fgType==FVM::FLUXGRIDTYPE_RADIAL)
                    return evaluateEnergyDistribution_fr(ir,p, dfdp, dfdr);
                else 
                    return evaluateEnergyDistribution(ir,p, dfdp, dfdr);
            }

            /**
             * Returns the non-normalized pitch distribution Fbar such that
             *  Fbar * F0
             * is proportional to the full distribution F up to some
             * arbitrary normalization factor.
             */
            virtual real_t evaluatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t *dfdxi0=nullptr, real_t *dfdp=nullptr, real_t *dfdr=nullptr) = 0;
            virtual real_t evaluatePitchDistribution_fr(len_t ir, real_t xi0, real_t p, real_t *dfdxi0=nullptr, real_t *dfdp=nullptr, real_t *dfdr=nullptr) = 0;
            virtual real_t evaluatePitchDistribution(len_t ir, real_t xi0, real_t p, FVM::fluxGridType fgType, real_t *dfdxi0=nullptr, real_t *dfdp=nullptr, real_t *dfdr=nullptr){
                if (fgType==FVM::FLUXGRIDTYPE_RADIAL)
                    return evaluatePitchDistribution_fr(ir,xi0,p,dfdxi0, dfdp, dfdr);
                else 
                    return evaluatePitchDistribution(ir,xi0,p,dfdxi0, dfdp, dfdr);   
            }

    };
}

#endif/*_DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_HPP*/