#ifndef _DREAM_FVM_NUMERIC_STELLARATOR_RADIAL_GRID_GENERATOR_HPP
#define _DREAM_FVM_NUMERIC_STELLARATOR_RADIAL_GRID_GENERATOR_HPP

#include <gsl/gsl_interp.h>
#include <softlib/SFile.h>
#include "FVM/Grid/Stellarator/RadialGridGeneratorStellarator.hpp"
#include "FVM/Grid/Stellarator/Interpolator3DSpatial.hpp"

namespace DREAM::FVM { 
    class NumericStellaratorRadialGridGenerator : public RadialGridGeneratorStellarator {
    public:
        struct eq_data {
            len_t nrho, ntheta, nphi;
            const real_t *rho, *theta, *phi;  // Coordinate arrays (1D)
            const real_t *dataG, *dataI;      // toroial and poloidal magnetic field strengths (1D)
            const real_t *dataiota;           // Rotational transform (1D)
            const real_t *datapsi;            // Toroidal flux (1D)

            const real_t *dataB, *dataBdotGradphi;           // Magnetic field data (3D)
            const real_t *dataJacobian, *datagtt, *datagtp;  // Jacobian data (3D)
            const real_t *datalambdat, *datalambdap;         // Stream function data (3D)
        };

    protected:

        len_t nphi_interp;
    
    private:
        struct eq_data *providedData;

        bool isBuilt = false;

        real_t rMin, rMax;
        real_t *rf_provided=nullptr;
        real_t *r=nullptr, *r_f=nullptr;

        // Input data
        len_t nrho, ntheta, nphi;

        // Interpolation objects for interpolating in input data
        gsl_spline *spline_G, *spline_I, *spline_iota,  *spline_psi;
        
        FVM::Interpolator3DSpatial *interp_B, *interp_Jacobian, //*interp_K, 
                            *interp_BdotGradphi, *interp_gtt, *interp_gtp,
                            *interp_lambdat, *interp_lambdap; 
        gsl_interp_accel *acc_r, *acc_theta;

		real_t *addThetaDataPoint(const real_t*, const len_t, const len_t, len_t);

        real_t _angleBounded(const real_t) const;

    public:
        NumericStellaratorRadialGridGenerator(
            const len_t nr, const real_t r0, const real_t ra, 
            const real_t R0, const len_t nfp, struct eq_data*,
			const len_t ntheta_interp=64, const len_t nphi_interp=64
        );
        NumericStellaratorRadialGridGenerator(
            const real_t *r_f, const len_t nr, 
            const real_t R0, const len_t nfp, struct eq_data*,
			const len_t ntheta_interp=64, const len_t nphi_interp=64
        );
        ~NumericStellaratorRadialGridGenerator();

        virtual bool NeedsRebuild(const real_t) const override { return (!isBuilt); }
        virtual bool Rebuild(const real_t, RadialGridStellarator*);

		real_t EvalB(const real_t, const real_t, const real_t);
        
        // Override virtual methods
		virtual real_t BAtThetaPhi(const len_t, const real_t, const real_t) override;
		virtual real_t BAtThetaPhi_f(const len_t, const real_t, const real_t) override;

        virtual real_t JacobianAtThetaPhi(const len_t ir, const real_t theta, const real_t phi) override { return JacobianAtThetaPhi(r[ir], theta, phi); }
        virtual real_t JacobianAtThetaPhi_f(const len_t ir, const real_t theta, const real_t phi) override { return JacobianAtThetaPhi(r_f[ir], theta, phi); }
        real_t JacobianAtThetaPhi(const real_t, const real_t, const real_t);

        virtual real_t BdotGradphiAtThetaPhi(const len_t ir, const real_t theta, const real_t phi) override { return BdotGradphiAtThetaPhi(r[ir], theta, phi); }
        virtual real_t BdotGradphiAtThetaPhi_f(const len_t ir, const real_t theta, const real_t phi) override { return BdotGradphiAtThetaPhi(r_f[ir], theta, phi); }
        real_t BdotGradphiAtThetaPhi(const real_t, const real_t, const real_t);

        virtual real_t gttAtThetaPhi(const len_t ir, const real_t theta, const real_t phi) override { return gttAtThetaPhi(r[ir], theta, phi); }
        virtual real_t gttAtThetaPhi_f(const len_t ir, const real_t theta, const real_t phi) override { return gttAtThetaPhi(r_f[ir], theta, phi); }
        real_t gttAtThetaPhi(const real_t, const real_t, const real_t);

        virtual real_t gtpAtThetaPhi(const len_t ir, const real_t theta, const real_t phi) override { return gtpAtThetaPhi(r[ir], theta, phi); }
        virtual real_t gtpAtThetaPhi_f(const len_t ir, const real_t theta, const real_t phi) override { return gtpAtThetaPhi(r_f[ir], theta, phi); }
        real_t gtpAtThetaPhi(const real_t, const real_t, const real_t);

        virtual void EvaluateGeometricQuantities(
            const len_t ir, const real_t theta, const real_t phi, real_t &B, real_t &Jacobian, real_t &BdotGradphi, real_t &gttOverJ2, real_t &gtpOverJ2
        )  override { EvaluateGeometricQuantities(r[ir], theta, phi, B, Jacobian, BdotGradphi, gttOverJ2, gtpOverJ2); }
        virtual void EvaluateGeometricQuantities_fr(
            const len_t ir, const real_t theta, const real_t phi, real_t &B, real_t &Jacobian, real_t &BdotGradphi, real_t &gttOverJ2, real_t &gtpOverJ2
        )  override { EvaluateGeometricQuantities(r_f[ir], theta, phi, B, Jacobian, BdotGradphi, gttOverJ2, gtpOverJ2); }

        void EvaluateGeometricQuantities(
            const real_t r, const real_t theta, real_t phi, real_t &B, real_t &Jacobian, real_t &BdotGradphi, real_t &gttOverJ2, real_t &gtpOverJ2
        );

		// Output generation helper routines
		virtual const len_t GetNPsi() override { return this->GetNr(); }
		virtual const len_t GetNTheta() override { return this->ntheta; }
        virtual const len_t GetNPhi() override { return this->nphi; }

		virtual const real_t *GetPoloidalAngle() override;
		virtual const real_t *GetToroidalAngle() override;

        // Debugging method
        void __SaveB(const char*); // TODO: Remove?
    };
}

#endif/*_DREAM_FVM_NUMERIC_STELLARATOR_RADIAL_GRID_GENERATOR_HPP*/ // TODO
