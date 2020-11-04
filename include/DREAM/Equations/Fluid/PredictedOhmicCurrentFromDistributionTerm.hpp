#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IonHandler.hpp"
/**
 * Implementation of a class which represents the sigma_numerical*E contribution 
 * that we expect f_hot to carry in steady-state for small E. 
 * 
 */
namespace DREAM {
    class PredictedOhmicCurrentFromDistributionTerm : public FVM::DiagonalComplexTerm {
    private:
        RunawayFluid *REFluid;
        IonHandler *ionHandler;
        real_t scaleFactor;
        len_t id_ni;
    protected:
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override {
            bool useCollisionlessLimit = true;
            len_t offset = 0;
            for(len_t n = 0; n<nMultiples; n++)
                for (len_t ir = 0; ir < nr; ir++){
                    real_t Zeff = ionHandler->evaluateZeff(ir);
                    real_t fracCond = FractionOfBraamsConductivity(Zeff);
                    real_t w0 = scaleFactor * fracCond / sqrt(grid->GetRadialGrid()->GetFSA_B2(ir));
                    real_t dSigmaBraams = REFluid->evaluatePartialContributionSauterConductivity(ir, derivId, n, useCollisionlessLimit);
                    for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                            diffWeights[offset + i] = w0*dSigmaBraams;
                    offset += n1[ir]*n2[ir];
                }
            // if ions, add the Zeff contribution from FractionOfBraams
            if(derivId == id_ni){
            len_t offset = 0;
                for(len_t n = 0; n<nMultiples; n++){
                    len_t iz,Z0;
                    ionHandler->GetIonIndices(n,iz,Z0);
                    for (len_t ir = 0; ir < nr; ir++){
                        real_t nfree = ionHandler->evaluateFreeElectronDensityFromQuasiNeutrality(ir);
                        real_t nZ0Z0 = ionHandler->evaluateZ0Z0(ir);
                        real_t Zeff = ionHandler->evaluateZeff(ir);
                        real_t dFracCond; // get the Zeff derivative of fracCond
                        if(nfree==0)
                            dFracCond = 0;
                        else {
                            FractionOfBraamsConductivity(Zeff,&dFracCond);
                            dFracCond *= Z0/nfree * (Z0 - nZ0Z0/nfree); // add dZeff/dni
                        }
                        real_t dw0 = scaleFactor * dFracCond / sqrt(grid->GetRadialGrid()->GetFSA_B2(ir));
                        real_t sigmaBraams = REFluid->evaluateSauterElectricConductivity(ir,useCollisionlessLimit);
                        for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                                diffWeights[offset + i] += dw0*sigmaBraams;
                        offset += n1[ir]*n2[ir];
                    }
                }
            }
        }

        // Set the weights
        virtual void SetWeights() override {
            len_t offset = 0;
            bool useCollisionlessLimit = true;
            for (len_t ir = 0; ir < nr; ir++){
                real_t Zeff = ionHandler->evaluateZeff(ir);
                real_t fracCond = FractionOfBraamsConductivity(Zeff);
                real_t sigma = REFluid->evaluateSauterElectricConductivity(ir, useCollisionlessLimit);
                real_t w = scaleFactor * fracCond / sqrt(grid->GetRadialGrid()->GetFSA_B2(ir));
                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                        weights[offset + i] = w * sigma;
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
        static real_t FractionOfBraamsConductivity(real_t Zeff, real_t *dFracDZ=nullptr){
            real_t a = -1.406;
            real_t b = 1.888;
            if(dFracDZ != nullptr) // output Zeff derivative of the return value
                *dFracDZ = -a/((b+Zeff)*(b+Zeff));

            return 1 + a/(b + Zeff);
        }
        PredictedOhmicCurrentFromDistributionTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u, RunawayFluid *ref, IonHandler *ih, real_t scaleFactor = 1.0) 
            : FVM::DiagonalComplexTerm(g,u), REFluid(ref), ionHandler(ih), scaleFactor(scaleFactor)
        {
            this->id_ni = u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

            /**
             * So far, we only account for the temperature dependence in the conductivity 
             * Jacobian and not, for example, ion densities which would enter through Zeff
             * and n_cold via the collisionality in the neoclassical corrections. 
             */
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD));
            AddUnknownForJacobian(unknowns,this->id_ni);
        }

    };
}
