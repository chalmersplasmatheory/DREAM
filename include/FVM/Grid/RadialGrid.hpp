#ifndef _DREAM_FVM_RADIAL_GRID_HPP
#define _DREAM_FVM_RADIAL_GRID_HPP

namespace DREAM::FVM { class RadialGrid; }

#include "FVM/FVMException.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include "FVM/Grid/FluxSurfaceAverager.hpp"
#include <functional> 

namespace DREAM::FVM {
	class RadialGrid {
    public:
        // Specification for functions used in flux surface averages
        static real_t FSA_FUNC_UNITY(real_t, real_t, real_t, void*)
            {return 1;};
        static real_t FSA_FUNC_ONE_OVER_R_SQUARED(real_t, real_t ROverR0,real_t, void*)
            {return 1/(ROverR0*ROverR0);}
        static real_t FSA_FUNC_B(real_t BOverBmin, real_t,real_t, void*)
            {return BOverBmin;}
        static real_t FSA_FUNC_B_SQUARED(real_t BOverBmin, real_t,real_t, void*)
            {return BOverBmin*BOverBmin;}
        static real_t FSA_FUNC_NABLA_R_SQUARED_OVER_R_SQUARED(real_t, real_t ROverR0,real_t NablaR2, void*)
            {return NablaR2/(ROverR0*ROverR0);}

        struct EPF_params {real_t x; real_t BminOverBmax; len_t ir; RadialGrid *rGrid; fluxGridType fgType;};
        static real_t FSA_FUNC_EFF_PASS_FRAC(real_t BOverBmin, real_t, real_t, void *par){ 
            struct EPF_params *params = (struct EPF_params *) par;
            real_t BminOverBmax = params->BminOverBmax;
            real_t x = params->x;
            return sqrt(1 - x * BminOverBmax * BOverBmin );
        }

        static real_t FSA_FUNC_XI(real_t BOverBmin, real_t, real_t, void *xiPtr){ 
            real_t xi0 = *(real_t*)xiPtr;
            if(BOverBmin < 1 + 100*realeps)
                return xi0;
            real_t xi = sqrt(1 - (1-xi0*xi0)*BOverBmin);
            if(xi0>=0)
                return xi;
            else 
                return -xi;
        }
        

        // Alternative parametric representation of FSA functions
        static constexpr int_t 
            FSA_PARAM_UNITY[4] = {0,0,0,1},
            FSA_PARAM_ONE_OVER_R_SQUARED[4] = {0,-2,0,1},
            FSA_PARAM_B[4] = {1,0,0,1},
            FSA_PARAM_B_SQUARED[4] = {2,0,0,1},
            FSA_PARAM_NABLA_R_SQUARED_OVER_R_SQUARED[4] = {0,-2,1,1};


        // Specification for functions used in bounce averages
        static real_t BA_FUNC_UNITY(real_t,real_t,real_t,real_t,void*)
            {return 1;}
        static real_t BA_FUNC_XI(real_t xiOverXi0,real_t,real_t,real_t,void*)
            {return xiOverXi0;}
        static real_t BA_FUNC_XI_SQUARED_OVER_B(real_t xiOverXi0,real_t BOverBmin,real_t,real_t,void*)
            {return xiOverXi0*xiOverXi0/BOverBmin;}
        static real_t BA_FUNC_B_CUBED(real_t, real_t BOverBmin, real_t, real_t, void*)
            {return BOverBmin*BOverBmin*BOverBmin;}
        static real_t BA_FUNC_XI_SQUARED_B_SQUARED(real_t xiOverXi0, real_t BOverBmin, real_t, real_t, void*)
            {return BOverBmin*BOverBmin*xiOverXi0*xiOverXi0;}

        // Alternative representation of functions to be bounce averaged:
        // lists containing exponents of the various contributing factors
        static constexpr int_t 
            BA_PARAM_UNITY[5] = {0,0,0,0,1},
            BA_PARAM_XI[5] = {1,0,0,0,1},
            BA_PARAM_XI_SQUARED_OVER_B[5] = {2,-1,0,0,1},
            BA_PARAM_B_CUBED[5] = {0,3,0,0,1},
            BA_PARAM_XI_SQUARED_B_SQUARED[5] = {2,2,0,0,1};

	private:
        // Flux-surface averaged quantities.
        real_t 
            *effectivePassingFraction   = nullptr, // Per's Eq (11.24)
            *effectivePassingFraction_f = nullptr, // Per's Eq (11.24)
            *FSA_B                      = nullptr, // <B> / Bmin
            *FSA_B_f                    = nullptr, // <B> / Bmin
            *FSA_B2                     = nullptr, // <B^2> / Bmin^2
            *FSA_B2_f                   = nullptr, // <B^2> / Bmin^2
            *FSA_nablaR2OverR2          = nullptr, // R0^2*<|nabla r|^2/R^2>
            *FSA_nablaR2OverR2_f        = nullptr, // R0^2*<|nabla r|^2/R^2>
            *FSA_1OverR2                = nullptr, // R0^2*<1/R^2>
            *FSA_1OverR2_f              = nullptr; // R0^2*<1/R^2>

        // Number of radial grid points
        len_t nr;

        // Radial grid
        // NOTE that 'r' has 'nr' elements, while
        // 'r_f' has (nr+1) elements.
        real_t *r=nullptr, *r_f=nullptr;
        // Radial grid steps
        //   dr[i]   = r_f[i+1] - r_f[i]   (nr elements)
        //   dr_f[i] = r[i+1] - r[i]       (nr-1 elements)
        real_t *dr=nullptr, *dr_f=nullptr;

        // Magnetic field quantities
        real_t 
            *Bmin       = nullptr,
            *Bmin_f     = nullptr,
            *Bmax       = nullptr,
            *Bmax_f     = nullptr,
            *xi0TrappedBoundary   = nullptr,
            *xi0TrappedBoundary_f = nullptr,
            *BtorGOverR0   = nullptr,
            *BtorGOverR0_f = nullptr,
            *psiPrimeRef   = nullptr,
            *psiPrimeRef_f = nullptr,
            *psiToroidal   = nullptr,
            *psiToroidal_f = nullptr,
            R0;
        
        // Orbit-phase-space Jacobian factors
        real_t
             *VpVol = nullptr,    // Size NR
             *VpVol_f = nullptr;  // Size NR+1

        // Deallocators
        void DeallocateReferenceMagneticData(){
            if(BtorGOverR0 == nullptr)
                return;
            delete [] BtorGOverR0;
            delete [] BtorGOverR0_f;
            delete [] psiPrimeRef;
            delete [] psiPrimeRef_f;
        }
        void DeallocateMagneticExtremumData(){
            if(Bmin == nullptr)
                return;
            delete [] Bmin;
            delete [] Bmin_f;
            delete [] Bmax;
            delete [] Bmax_f;
            delete [] xi0TrappedBoundary;
            delete [] xi0TrappedBoundary_f;
        }
        void SetFluxSurfaceAverage(real_t *&FSA_quantity, real_t *&FSA_quantity_f, real_t(*F)(real_t,real_t,real_t,void*), void *par=nullptr, const int_t *Flist = nullptr);

        virtual void RebuildFluxSurfaceAveragedQuantities();
        void SetEffectivePassingFraction(real_t*&, real_t*&, real_t*, real_t*);
        static real_t effectivePassingFractionIntegrand(real_t x, void *p);

        void DeallocateFSAvg();
        void InitializeFSAvg(
            real_t *epf, real_t *epf_f, real_t *Bavg, real_t *Bavg_f, 
            real_t *B2avg, real_t *B2avg_f,
            real_t *OneOverR2_avg, real_t *OneOverR2_avg_f,
            real_t *nablaR2OverR2_avg, real_t *nablaR2OverR2_avg_f);

        static constexpr real_t realeps = std::numeric_limits<real_t>::epsilon();    

	protected:
        FluxSurfaceAverager *fluxSurfaceAverager;
        RadialGridGenerator *generator;

    public:
        RadialGrid(RadialGridGenerator*, const real_t t0=0, 
            FluxSurfaceAverager::interp_method im = FluxSurfaceAverager::INTERP_STEFFEN,
            FluxSurfaceAverager::quadrature_method qm = FluxSurfaceAverager::QUAD_FIXED_LEGENDRE
        );
        virtual ~RadialGrid();

        void DeallocateGrid();

        void Initialize(
            real_t *r, real_t *r_f,
            real_t *dr, real_t *dr_f
        ) {
            if (this->r == r &&
                this->r_f == r_f &&
                this->dr == dr &&
                this->dr_f == dr_f)
                return;

            DeallocateGrid();

            this->r    = r;
            this->r_f  = r_f;
            this->dr   = dr;
            this->dr_f = dr_f;
        }

        void SetReferenceMagneticFieldData(
            real_t *BtorGOverR0, real_t *BtorGOverR0_f,
            real_t *psiPrimeRef, real_t *psiPrimeRef_f,
            real_t R0
        );
        void SetMagneticExtremumData(
            real_t *Bmin, real_t *Bmin_f,
            real_t *Bmax, real_t *Bmax_f,
            real_t *theta_Bmin, real_t *theta_Bmin_f,
            real_t *theta_Bmax, real_t *theta_Bmax_f,
            real_t *xi0TrappedBoundary, real_t *xi0TrappedBoundary_f
        );

        bool Rebuild(const real_t);

        virtual void RebuildJacobians();
        
        real_t CalculateFluxSurfaceAverage(len_t ir, fluxGridType fluxGridType, real_t(*F)(real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);
        real_t EvaluateFluxSurfaceIntegral(len_t ir, fluxGridType fluxGridType, real_t(*F)(real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);
        real_t CalculatePXiBounceAverageAtP(len_t ir, real_t xi0, fluxGridType fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);
        real_t EvaluatePXiBounceIntegralAtP(len_t ir, real_t xi0, fluxGridType fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);
        void SetVpVol(real_t *VpVol, real_t *VpVol_f){
            if(this->VpVol!=nullptr){
                delete [] this->VpVol;
                delete [] this->VpVol_f;
            }        
            this->VpVol = VpVol;
            this->VpVol_f = VpVol_f;
        }

        void GetRThetaPhiFromCartesian(real_t *r, real_t *theta, real_t *phi, real_t x, real_t y, real_t z, real_t lengthScale, real_t startingGuessR){return this->generator->GetRThetaPhiFromCartesian(r,theta,phi,x,y,z,lengthScale, startingGuessR);}
        void GetGradRCartesian(real_t *gradRCartesian, real_t r, real_t theta, real_t phi){return this->generator->GetGradRCartesian(gradRCartesian,r,theta, phi);}
        real_t FindClosestApproach(real_t x1, real_t y1, real_t z1, real_t x2, real_t y2, real_t z2){
            return this->generator->FindClosestApproach(x1, y1, z1, x2, y2, z2);
        }

        /**
         * Getters of magnetic field strength quantities
         */
        const real_t *GetBmin() const {return this->Bmin;}
        const real_t  GetBmin(const len_t ir) const {return this->Bmin[ir];}
        const real_t *GetBmin_f() const {return this->Bmin_f;}
        const real_t  GetBmin_f(const len_t ir) const {return this->Bmin_f[ir];}
        const real_t *GetBmax() const {return this->Bmax;}
        const real_t  GetBmax(const len_t ir) const {return this->Bmax[ir];}
        const real_t *GetBmax_f() const {return this->Bmax_f;}
        const real_t  GetBmax_f(const len_t ir) const {return this->Bmax_f[ir];}
        const real_t *GetBTorG() const {return this->BtorGOverR0;}
        const real_t  GetBTorG(const len_t ir) const {return this->BtorGOverR0[ir];}
        const real_t *GetBTorG_f() const {return this->BtorGOverR0_f;}
        const real_t  GetBTorG_f(const len_t ir) const {return this->BtorGOverR0_f[ir];}
        
        // Returns the xi0 value corresponding to the positive 
        // trapped-passing boundary at radial index ir
        const real_t GetXi0TrappedBoundary(const len_t ir) const 
            {return xi0TrappedBoundary[ir];}
        const real_t* GetXi0TrappedBoundary() const 
            {return xi0TrappedBoundary;}
        // Returns trapped-passing boundary on radial flux grid
        const real_t GetXi0TrappedBoundary_fr(const len_t ir) const 
            {return xi0TrappedBoundary_f[ir];}
        const real_t* GetXi0TrappedBoundary_fr() const 
            {return xi0TrappedBoundary_f;}

        /**
         * Getters of grid data:
         */
        // Returns the number of radial grid points in this grid
        len_t GetNr() const { return this->nr; }
        // Returns the vector containing all radial grid points
        const real_t *GetR() const { return this->r; }
        const real_t  GetR(const len_t i) const { return this->r[i]; }
        const real_t *GetR_f() const { return this->r_f; }
        const real_t  GetR_f(const len_t i) const { return this->r_f[i]; }
        const real_t  GetR0() const { return this->R0;}
        // Returns a vector containing all radial steps
        const real_t *GetDr() const { return this->dr; }
        const real_t  GetDr(const len_t i) const { return this->dr[i]; }
        const real_t *GetDr_f() const { return this->dr_f; }
        const real_t  GetDr_f(const len_t i) const { return this->dr_f[i]; }
        
        const real_t GetMinorRadius() const { return r_f[this->nr]; }

		// Routines used for saving equilibrium to output file
		virtual const real_t GetZ0() { return this->generator->GetZ0(); }
		virtual const len_t GetNPsi() { return this->generator->GetNPsi(); }
		virtual const len_t GetNTheta() { return this->generator->GetNTheta(); }
		virtual const real_t *GetFluxSurfaceRMinusR0() { return this->generator->GetFluxSurfaceRMinusR0(); }
		virtual const real_t *GetFluxSurfaceRMinusR0_f() { return this->generator->GetFluxSurfaceRMinusR0_f(); }
		virtual const real_t *GetFluxSurfaceZMinusZ0() { return this->generator->GetFluxSurfaceZMinusZ0(); }
		virtual const real_t *GetFluxSurfaceZMinusZ0_f() { return this->generator->GetFluxSurfaceZMinusZ0_f(); }
		virtual const real_t *GetPoloidalAngle() { return this->generator->GetPoloidalAngle(); }
        
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
        const real_t SafetyFactorNormalized(const len_t ir, const real_t mu0Ip) const {
            if(mu0Ip==0)
                return std::numeric_limits<real_t>::infinity();
            real_t twoPi = 2*M_PI;
            real_t twoPiCubed = twoPi*twoPi*twoPi;
            return VpVol[ir]*VpVol[ir]/(twoPiCubed*mu0Ip) * GetBTorG(ir)  
                    * GetFSA_1OverR2(ir) * GetFSA_NablaR2OverR2(ir);
        }

        const real_t *GetToroidalFlux() const 
            { return psiToroidal; }
        const real_t GetToroidalFlux(len_t ir) const 
            { return psiToroidal[ir]; }
        const real_t *GetToroidalFlux_f() const
            { return psiToroidal_f; }
        const real_t GetToroidalFlux_f(len_t ir) const
            { return psiToroidal_f[ir]; }
        
        /**
         * Getters of flux-surface averaged Jacobian
         */
        const real_t *GetVpVol() const {return this->VpVol; }
        const real_t  GetVpVol(const len_t ir) const {return this->VpVol[ir]; }
        const real_t *GetVpVol_f() const {return this->VpVol_f; }
        const real_t  GetVpVol_f(const len_t ir) const {return this->VpVol_f[ir]; }
        
        /**
         * Getters of flux surface averaged quantities
         */
        const real_t  *GetEffPassFrac() const { return this->effectivePassingFraction; }
        const real_t   GetEffPassFrac(const len_t ir) const { return this->effectivePassingFraction[ir]; }
        const real_t  *GetFSA_B2() const { return this->FSA_B2; }
        const real_t   GetFSA_B2(const len_t ir) const { return this->FSA_B2[ir]; }
        const real_t  *GetFSA_B2_f() const { return this->FSA_B2_f; }
        const real_t   GetFSA_B2_f(const len_t ir) const { return this->FSA_B2_f[ir]; }
        const real_t  *GetFSA_B() const { return this->FSA_B; }
        const real_t   GetFSA_B(const len_t ir) const { return this->FSA_B[ir]; }
        const real_t  *GetFSA_B_f() const { return this->FSA_B_f; }
        const real_t   GetFSA_B_f(const len_t ir) const { return this->FSA_B_f[ir]; }
        const real_t  *GetFSA_1OverR2() const { return this->FSA_1OverR2; }
        const real_t   GetFSA_1OverR2(const len_t ir) const { return this->FSA_1OverR2[ir]; }
        const real_t  *GetFSA_1OverR2_f() const { return this->FSA_1OverR2_f; }
        const real_t   GetFSA_1OverR2_f(const len_t ir) const { return this->FSA_1OverR2_f[ir]; }
        const real_t  *GetFSA_NablaR2OverR2() const { return this->FSA_nablaR2OverR2; }
        const real_t   GetFSA_NablaR2OverR2(const len_t ir) const { return this->FSA_nablaR2OverR2[ir]; }
        const real_t  *GetFSA_NablaR2OverR2_f() const { return this->FSA_nablaR2OverR2_f; }
        const real_t   GetFSA_NablaR2OverR2_f(const len_t ir) const { return this->FSA_nablaR2OverR2_f[ir]; }

        FluxSurfaceAverager *GetFluxSurfaceAverager(){return fluxSurfaceAverager;}

        bool NeedsRebuild(const real_t t) const { return this->generator->NeedsRebuild(t); }

	};

    class RadialGridException : public FVMException {
    public:
        template<typename ... Args>
        RadialGridException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("RadialGrid");
        }
    };
}

#endif/*_DREAM_FVM_RADIAL_GRID_HPP*/
