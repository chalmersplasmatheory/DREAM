#ifndef _DREAM_FVM_NUMERIC_B_RADIAL_GRID_GENERATOR_HPP // TODO
#define _DREAM_FVM_NUMERIC_B_RADIAL_GRID_GENERATOR_HPP // TODO

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline2d.h>
#include <string>
#include <softlib/SFile.h>
#include "FVM/Grid/RadialGridGenerator.hpp"
#include "FVM/Grid/Stellarator/Interpolator3D.hpp"

namespace DREAM::FVM {
    class NumericBRadialGridGenerator : public RadialGridGenerator {
    public:
        enum file_format {
            FILE_FORMAT_DESC
        };

    protected:
        virtual real_t FindMagneticFieldExtremum(
			len_t ir, int_t sgn, enum fluxGridType
		) override;
        gsl_multimin_fminimizer *gsl_multi_fmin;

        len_t nphi_interp;
        real_t  *BpolIOverR0,        *BpolIOverR0_f;        // poloidal magnetic field strength (I(r)/R0)
        real_t  *iota,               *iota_f;               // Rotational transform
    
    private:
        bool isBuilt = false;

        real_t rMin, rMax;
        real_t *rf_provided=nullptr;
        real_t *r=nullptr, *r_f=nullptr;

        // Input data
        real_t R0;
        len_t npsi, ntheta, nphi;
        real_t *input_r=nullptr;    // Calculated minor radius in outer midplane (same size as psi)
        real_t *psi=nullptr, *theta=nullptr, *phi=nullptr;        // Poloidal flux, poloidal angle and toroidal angle grids (1D)
        real_t *dataG=nullptr, *dataI=nullptr;// Magnetic field components (1D)
        real_t *dataiota=nullptr; // Rotational transform (1D)
        real_t //*dataK=nullptr
                *dataBdotGradphi=nullptr, *datagtt=nullptr, *datagtp=nullptr;     // Magnetic field and Jacobian data (3D)
        real_t *dataB=nullptr, *dataJacobian=nullptr;   // Magnetic field and Jacobian data (2D)

        std::string name;

        // Interpolation objects for interpolating in input data
        gsl_spline *spline_psi, *spline_G, *spline_I, *spline_iota;
        gsl_spline2d;
            //*spline_B, *spline_Jacobian;
        
        FVM::Interpolator3D *interp_B, *interp_Jacobian, //*interp_K, 
                            *interp_BdotGradphi, *interp_gtt, *interp_gtp; 
        gsl_interp_accel *acc_r, *acc_theta;

		real_t *addThetaDataPoint(const real_t*, const len_t, const len_t, len_t);

        void Init(const std::string&, enum file_format, const len_t, const len_t);

        real_t _angleBounded(const real_t) const;

    public:
        NumericBRadialGridGenerator(
            const len_t nr, const real_t r0, const real_t ra,
            const std::string&, enum file_format frmt=FILE_FORMAT_DESC,
			const len_t ntheta_interp=30, const len_t nphi_interp=30 // TODO: What is a good resolution?
        );
        NumericBRadialGridGenerator(
            const real_t *r_f, const len_t nr, const std::string&,
            enum file_format frmt=FILE_FORMAT_DESC,
			const len_t ntheta_interp=30, const len_t nphi_interp=30 // TODO: What is a good resolution?
        );
        ~NumericBRadialGridGenerator();

        len_t GetNphiInterp() const { return this->nphi_interp; }

        void LoadMagneticFieldData(const std::string&, enum file_format frmt=FILE_FORMAT_DESC);
        void LoadMagneticFieldData(SFile*, enum file_format frmt=FILE_FORMAT_DESC);

        virtual bool NeedsRebuild(const real_t) const override { return (!isBuilt); }
        virtual bool Rebuild(const real_t, RadialGrid*) override;
        /* TODO: If B has several min and max on flux-surface, we have to redo this 
        virtual void RebuildJacobians(RadialGrid*) override;*/

		real_t EvalB(const real_t, const real_t);
        
        // Override virtual methods
		virtual real_t BAtThetaPhi(const len_t, const real_t, const real_t) override;
		virtual real_t BAtThetaPhi_f(const len_t, const real_t, const real_t) override;

        virtual real_t JacobianAtThetaPhi(const len_t ir, const real_t theta, const real_t phi) { return JacobianAtThetaPhi(r[ir], theta, phi); }
        virtual real_t JacobianAtThetaPhi_f(const len_t ir, const real_t theta, const real_t phi) { return JacobianAtThetaPhi(r_f[ir], theta, phi); };
        real_t JacobianAtThetaPhi(const real_t, const real_t, const real_t);

        virtual real_t BdotGradphiAtThetaPhi(const len_t ir, const real_t theta, const real_t phi) { return BdotGradphiAtThetaPhi(r[ir], theta, phi); }
        virtual real_t BdotGradphiAtThetaPhi_f(const len_t ir, const real_t theta, const real_t phi) { return BdotGradphiAtThetaPhi(r_f[ir], theta, phi); };
        real_t BdotGradphiAtThetaPhi(const real_t, const real_t, const real_t);

        virtual real_t gttAtThetaPhi(const len_t ir, const real_t theta, const real_t phi) { return gttAtThetaPhi(r[ir], theta, phi); }
        virtual real_t gttAtThetaPhi_f(const len_t ir, const real_t theta, const real_t phi) { return gttAtThetaPhi(r_f[ir], theta, phi); };
        real_t gttAtThetaPhi(const real_t, const real_t, const real_t);

        virtual real_t gtpAtThetaPhi(const len_t ir, const real_t theta, const real_t phi) { return gtpAtThetaPhi(r[ir], theta, phi); }
        virtual real_t gtpAtThetaPhi_f(const len_t ir, const real_t theta, const real_t phi) { return gtpAtThetaPhi(r_f[ir], theta, phi); };
        real_t gtpAtThetaPhi(const real_t, const real_t, const real_t);

        virtual void EvaluateGeometricQuantities(
            const len_t ir, const real_t theta, const real_t phi, real_t &B, real_t &Jacobian, real_t &BdotGradphi, real_t &gttOverJ2, real_t &gtpOverJ2
        ) override { EvaluateGeometricQuantities(r[ir], theta, phi, B, Jacobian, BdotGradphi, gttOverJ2, gtpOverJ2); }
        virtual void EvaluateGeometricQuantities_fr(
            const len_t ir, const real_t theta, const real_t phi, real_t &B, real_t &Jacobian, real_t &BdotGradphi, real_t &gttOverJ2, real_t &gtpOverJ2
        ) override { EvaluateGeometricQuantities(r_f[ir], theta, phi, B, Jacobian, BdotGradphi, gttOverJ2, gtpOverJ2); }

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
        void __SaveB(const char*);
    };

    // TODO: Make DESC relevant
    struct NumericBData {
        std::string name;
        sfilesize_t npsi, ntheta, nphi;
        double R0;
        double *psi=nullptr, *theta=nullptr, *phi=nullptr;    // Coordinate grids
        double *G=nullptr, *I=nullptr, iota=nullptr; // Dim npsi
        double //*K=nullptr, 
                *BdotGradphi=nullptr; // Dim: npsi*ntheta*nphi? Is K even needed?
        double *B=nullptr; // Dim npsi*ntheta? Or npsi*ntheta*nphi? Should be independent of phi
        double *Jacobian=nullptr; // Dim npsi*ntheta? Or npsi*ntheta*nphi? Should be independent of phi
        double *gtt=nullptr, *gtp=nullptr; // Dim npsi*ntheta*nphi?

        ~NumericBData() {
            if (psi!=nullptr) delete [] psi;
            if (theta!=nullptr) delete [] theta;
            if (phi!=nullptr) delete [] phi;
            if (G!=nullptr) delete [] G;
            if (I!=nullptr) delete [] I;
            if (iota!=nullptr) delete [] iota;
            //if (K!=nullptr) delete [] K;
            if (BdotGradphi!=nullptr) delete [] BdotGradphi;
            if (B!=nullptr) delete [] B;
            if (Jacobian!=nullptr) delete [] Jacobian;
            if (gtt!=nullptr) delete [] gtt;
            if (gtp!=nullptr) delete [] gtp;
        }
    };

    struct NumericBData *LoadNumericBFromDESC(SFile*, const std::string& path="");
}

#endif/*_DREAM_FVM_NUMERIC_B_RADIAL_GRID_GENERATOR_HPP*/ // TODO
