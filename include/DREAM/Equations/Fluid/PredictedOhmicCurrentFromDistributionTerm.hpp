#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IonHandler.hpp"
/**
 * Implementation of a class which represents the sigma_numerical*E contribution 
 * that we expect f_hot to carry in steady-state for small E. 
 * TODO: add trapping correction
 */
namespace DREAM {
    class PredictedOhmicCurrentFromDistributionTerm : public FVM::DiagonalComplexTerm {
    private:
        RunawayFluid *REFluid;
        IonHandler *ionHandler;
        real_t scaleFactor;
    protected:
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override {
            real_t *dSigma = REFluid->evaluatePartialContributionBraamsConductivity(ionHandler->evaluateZeff(),derivId);

            len_t offset = 0;
            for(len_t n = 0; n<nMultiples; n++){
                for (len_t ir = 0; ir < nr; ir++){
                    real_t Zeff = ionHandler->evaluateZeff(ir);
                    real_t neoclassicalCorrection = REFluid->evaluateNeoclassicalConductivityCorrection(ir, Zeff, true);
                    real_t fracCond = neoclassicalCorrection * FractionOfBraamsConductivity(Zeff);
                    real_t w0 = scaleFactor * fracCond * sqrt(grid->GetRadialGrid()->GetFSA_B2(ir));
                    for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                            diffWeights[offset + i] = w0*dSigma[offset + i];
                    offset += n1[ir]*n2[ir];
                }
            }
        }

        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                //real_t w=0;
                real_t Zeff = ionHandler->evaluateZeff(ir);
                real_t neoclassicalCorrection = REFluid->evaluateNeoclassicalConductivityCorrection(ir, Zeff, true);
                real_t fracCond = neoclassicalCorrection * FractionOfBraamsConductivity(Zeff);
                real_t w = scaleFactor * fracCond * sqrt(grid->GetRadialGrid()->GetFSA_B2(ir)) * REFluid->evaluateBraamsElectricConductivity(ir,Zeff);
                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                        weights[offset + i] = w;
                offset += n1[ir]*n2[ir];
            }
        }
    public:
        
        /** 
         * A function of the form Frac = 1 + a / (b+Zeff) was
         * fitted to calculated conductivity ratios with CODE,
         * sigma_CODE / sigma_Braams, producing the below 
         * values for the parameters a and b. The 95% confidence
         * intervals give less than 5% variation in a and 10%
         * in b, and the maximum relative error in the fit was
         * 1%. 
         */
        static real_t FractionOfBraamsConductivity(real_t Zeff){
            real_t a = -1.406;
            real_t b = 1.888;

            return 1 + a/(b + Zeff);
        }
        PredictedOhmicCurrentFromDistributionTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u, RunawayFluid *ref, IonHandler *ih, real_t scaleFactor = 1.0) 
            : FVM::DiagonalComplexTerm(g,u), REFluid(ref), ionHandler(ih), scaleFactor(scaleFactor)
        {
            /**
             * So far, we only account for the temperature dependence in the conductivity 
             * Jacobian and not, for example, ion densities which would enter through Zeff
             * and n_cold via the collisionality in the neoclassical corrections. 
             */
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD));
        }

    };
}