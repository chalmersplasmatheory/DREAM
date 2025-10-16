#ifndef _DREAM_FVM_RADIAL_GRID_STELLARATOR_HPP // TODO
#define _DREAM_FVM_RADIAL_GRID_STELLARATOR_HPP

namespace DREAM::FVM { class RadialGridStellarator; }

#include "FVM/FVMException.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/Stellarator/RadialGridGeneratorStellarator.hpp"
#include "FVM/Grid/Stellarator/FluxSurfaceAveragerStellarator.hpp"
#include <functional>

namespace DREAM::FVM {
	class RadialGridStellarator : public RadialGrid {
    public:
        /** TODO: Take back code for BA, see original RadialGrid */


        // Specification for functions used in flux surface averages
        static real_t FSA_FUNC_UNITY(real_t, real_t, real_t, real_t, void*)
            {return 1;};
        static real_t FSA_FUNC_B(real_t BOverBmin, real_t, real_t, real_t, void*)
            {return BOverBmin;}
        static real_t FSA_FUNC_B_SQUARED(real_t BOverBmin, real_t, real_t, real_t, void*)
            {return BOverBmin*BOverBmin;}
        static real_t FSA_FUNC_1OverB(real_t BOverBmin, real_t, real_t, real_t, void*)
            {return 1 / BOverBmin;}
        static real_t FSA_FUNC_B_DOT_GRAD_PHI(real_t, real_t BdotGradphi, real_t, real_t, void*)
            {return BdotGradphi;}
        static real_t FSA_FUNC_GTT_OVER_JACOBIAN_SQUARED(real_t, real_t, real_t gttOverJ2, real_t, void*)
            {return gttOverJ2;}
        static real_t FSA_FUNC_GTP_OVER_JACOBIAN_SQUARED(real_t, real_t, real_t, real_t gtpOverJ2, void*)
            {return gtpOverJ2;}

        struct EPF_params {real_t x; real_t BminOverBmax; len_t ir; RadialGridStellarator *rGrid; fluxGridType fgType;};
        static real_t FSA_FUNC_EFF_PASS_FRAC(real_t BOverBmin, real_t, real_t, real_t, void *par){
            struct EPF_params *params = (struct EPF_params *) par;
            real_t BminOverBmax = params->BminOverBmax;
            real_t x = params->x;
            return sqrt(1 - x * BminOverBmax * BOverBmin );
        }


        // Alternative parametric representation of FSA functions
        static constexpr int_t
            FSA_PARAM_UNITY[5] = {0,0,0,0,1},
            FSA_PARAM_B[5] = {1,0,0,0,1},
            FSA_PARAM_B_SQUARED[5] = {2,0,0,0,1},
            FSA_PARAM_1OverB[5] = {-1,0,0,0,1},
            FSA_PARAM_B_DOT_GRAD_PHI[5] = {0,1,0,0,1},
            FSA_PARAM_GTT_OVER_JACOBIAN_SQUARED[5] = {0,0,1,0,1},
            FSA_PARAM_GTP_OVER_JACOBIAN_SQUARED[5] = {0,0,0,1,1};

        
	private:
        // Flux-surface averaged quantities.
        real_t
            *FSA_BdotGradphi            = nullptr, // <|B\dot\nabla\varphi|>
            *FSA_BdotGradphi_f          = nullptr, // <|B\dot\nabla\varphi|>
            *FSA_gttOverJ2              = nullptr, // <g_{\theta\theta}/J^2>
            *FSA_gttOverJ2_f            = nullptr, // <g_{\theta\theta}/J^2>
            *FSA_gtpOverJ2              = nullptr, // <g_{\theta\varphi}/J^2>
            *FSA_gtpOverJ2_f            = nullptr; // <g_{\theta\varphi}/J^2>

        // Magnetic field quantities
        real_t
            *BpolIOverR0   = nullptr,
            *BpolIOverR0_f = nullptr,
            *iota   = nullptr,
            *iota_f = nullptr,
            R0;

        // Orbit-phase-space Jacobian factors
        real_t
             *VpVol = nullptr,    // Size NR
             *VpVol_f = nullptr;  // Size NR+1

        // Deallocators
        void DeallocateReferenceMagneticData(){
            if(BpolIOverR0 == nullptr)
                return;
            delete [] BpolIOverR0;
            delete [] BpolIOverR0_f;
        }
        void DeallocateStellaratorData(){
            if(iota == nullptr)
                return;
            delete [] iota;
            delete [] iota_f;
        }
        void SetFluxSurfaceAverage(real_t *&FSA_quantity, real_t *&FSA_quantity_f, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par=nullptr, const int_t *Flist = nullptr);

        virtual void RebuildFluxSurfaceAveragedQuantities() override;
        void SetEffectivePassingFraction(real_t*&, real_t*&, real_t*, real_t*);
        static real_t effectivePassingFractionIntegrand(real_t x, void *p);

        void DeallocateFSAvg();
        void InitializeFSAvg(
            real_t *epf, real_t *epf_f, real_t *Bavg, real_t *Bavg_f, 
            real_t *B2avg, real_t *B2avg_f,
            real_t *OneOverB_avg, real_t *OneOverB_avg_f,
            real_t *BdotGradphi_avg, real_t *BdotGradphi_avg_f,
            real_t *gttOverJ2_avg, real_t *gttOverJ2_avg_f,
            real_t *gtpOverJ2_avg, real_t *gtpOverJ2_avg_f
        );

        static constexpr real_t realeps = std::numeric_limits<real_t>::epsilon();

	protected:
        FluxSurfaceAveragerStellarator *fluxSurfaceAverager;
        RadialGridGeneratorStellarator *generator;

    public:
        RadialGridStellarator(RadialGridGeneratorStellarator*, const real_t t0=0,
            FluxSurfaceAveragerStellarator::interp_method im = FluxSurfaceAveragerStellarator::INTERP_STEFFEN,
            FluxSurfaceAveragerStellarator::quadrature_method qm = FluxSurfaceAveragerStellarator::QUAD_FIXED_LEGENDRE
        );
        virtual ~RadialGridStellarator();

        void SetReferenceMagneticFieldData(
            real_t *BtorGOverR0, real_t *BtorGOverR0_f,
            real_t *BpolIOverR0, real_t *BpolIOverR0_f,
            real_t *psiPrimeRef, real_t *psiPrimeRef_f,
            real_t R0
        );
        void SetStellaratorData(
            real_t *iota, real_t *iota_f
        );
        void SetMagneticExtremumData(
            real_t *Bmin, real_t *Bmin_f,
            real_t *Bmax, real_t *Bmax_f,
            real_t *theta_Bmin, real_t *theta_Bmin_f,
            real_t *theta_Bmax, real_t *theta_Bmax_f,
            real_t *xi0TrappedBoundary, real_t *xi0TrappedBoundary_f
        );

        bool Rebuild(const real_t);

        virtual void RebuildJacobians() override;

        real_t CalculateFluxSurfaceAverage(len_t ir, fluxGridType fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);
        real_t EvaluateFluxSurfaceIntegral(len_t ir, fluxGridType fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);

        virtual const real_t *GetBPolI() const override {return this->BpolIOverR0;}
        virtual const real_t  GetBPolI(const len_t ir) const override {return this->BpolIOverR0[ir];}
        virtual const real_t *GetBPolI_f() const override {return this->BpolIOverR0_f;}
        virtual const real_t  GetBPolI_f(const len_t ir) const override {return this->BpolIOverR0_f[ir];}
		virtual const real_t *GetIota() const override {return this->iota;}
		virtual const real_t  GetIota(const len_t ir) const override {return this->iota[ir];}
		virtual const real_t *GetIota_f() const override {return this->iota_f;}
		virtual const real_t  GetIota_f(const len_t ir) const override {return this->iota_f[ir];}

		// Routines used for saving equilibrium to output file
        // TODO: Don't save equilibrium in output, or save something else
        virtual const real_t *GetToroidalAngle() { return this->generator->GetToroidalAngle(); }
        
        /**
         * Returns q*R0 on the distribution grid where q
         * is the safety factor and R0 the major radius.
         * The safety factor is signed, so that negative
         * currents yield negative safety factor, which keeps
         * track of the handed-ness of the field line twist
         *  ir: radial grid index
         *  mu0Ip: product of vacuum permeability and toroidal plasma
         *         current enclosed by the flux surface ir.
         */
        const real_t RotationalTransform(const len_t ir, const real_t mu0Ip) const {
            real_t twoPi = 2*M_PI;
            real_t twoPiCubed = twoPi*twoPi*twoPi;
            return twoPiCubed / (VpVol[ir]*VpVol[ir]*R0*R0)
                        * (mu0Ip / GetFSA_BdotGradphi(ir) - GetFSA_gtpOverJ2(ir)) 
                            / GetFSA_gttOverJ2(ir);
        }

        const real_t SafetyFactorNormalized(const len_t ir, const real_t /*mu0Ip*/) const {
            real_t iota = this->iota[ir]; 
            // TODO, future option, 
            // real_t iota = RotationalTransform(ir, mu0Ip); 
            return (this->BtorGOverR0[ir] + iota * BpolIOverR0[ir]) * R0 / iota * FSA_1OverB[ir] / Bmin[ir];
        }

        /**
         * Getters of flux surface averaged quantities
         */
        virtual const real_t  *GetFSA_BdotGradphi() const override { return this->FSA_BdotGradphi; }
        virtual const real_t   GetFSA_BdotGradphi(const len_t ir) const override { return this->FSA_BdotGradphi[ir]; }
        virtual const real_t  *GetFSA_BdotGradphi_f() const override { return this->FSA_BdotGradphi_f; }
        virtual const real_t   GetFSA_BdotGradphi_f(const len_t ir) const override { return this->FSA_BdotGradphi_f[ir]; }
        virtual const real_t  *GetFSA_gttOverJ2() const override { return this->FSA_gttOverJ2; }
        virtual const real_t   GetFSA_gttOverJ2(const len_t ir) const override { return this->FSA_gttOverJ2[ir]; }
        virtual const real_t  *GetFSA_gttOverJ2_f() const override { return this->FSA_gttOverJ2_f; }
        virtual const real_t   GetFSA_gttOverJ2_f(const len_t ir) const override { return this->FSA_gttOverJ2_f[ir]; }
        virtual const real_t  *GetFSA_gtpOverJ2() const override { return this->FSA_gtpOverJ2; }
        virtual const real_t   GetFSA_gtpOverJ2(const len_t ir) const override { return this->FSA_gtpOverJ2[ir]; }
        virtual const real_t  *GetFSA_gtpOverJ2_f() const override { return this->FSA_gtpOverJ2_f; }
        virtual const real_t   GetFSA_gtpOverJ2_f(const len_t ir) const override { return this->FSA_gtpOverJ2_f[ir]; }

        FluxSurfaceAveragerStellarator *GetFluxSurfaceAverager(){return fluxSurfaceAverager;}

	};

    class RadialGridStellaratorException : public FVMException {
    public:
        template<typename ... Args>
        RadialGridStellaratorException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("RadialGridStellarator");
        }
    };
}

#endif/*_DREAM_FVM_RADIAL_GRID_STELLARATOR_HPP*/
