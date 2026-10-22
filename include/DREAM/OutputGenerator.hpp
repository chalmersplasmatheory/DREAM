#ifndef _DREAM_OUTPUT_GENERATOR_HPP
#define _DREAM_OUTPUT_GENERATOR_HPP

namespace DREAM { class OutputGenerator; }

#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/OtherQuantityHandler.hpp"
#include "DREAM/Equations/Kinetic/MollerDeltaAngleKernel.hpp"

namespace DREAM {
	class OutputGenerator {
	public:
		struct static_qtys {
			MollerKernelHandler *mollerKernelHandler=nullptr;
		};

	protected:
		FVM::Grid *scalarGrid, *fluidGrid, *hottailGrid, *runawayGrid;
		FVM::UnknownQuantityHandler *unknowns;
		IonHandler *ions;
		OtherQuantityHandler *oqty;
		EquationSystem *eqsys;
		struct static_qtys *staticQuantities;

        bool savesettings = true;

		virtual void SaveCodeInfo(const std::string&) = 0;
		virtual void SaveGrids(const std::string&, bool) = 0;
		virtual void SaveGridBoltzmann(const std::string&) = 0;
		virtual void SaveGridBoltzmannDeltaKernel(const std::string&, const MollerDeltaAngleKernel*) = 0;
		virtual void SaveIonMetaData(const std::string&) = 0;
		virtual void SaveOtherQuantities(const std::string&) = 0;
		virtual void SaveSettings(const std::string&) = 0;
        virtual void SaveSolverData(const std::string&) = 0;
		virtual void SaveTimings(const std::string&) = 0;
		virtual void SaveUnknowns(const std::string&, bool) = 0;
	public:
		OutputGenerator(
			EquationSystem*, struct static_qtys *staticQuantities=nullptr,
			bool savesettings=true
		);
        virtual ~OutputGenerator();

		virtual void Save(bool current=false);
        virtual void SaveCurrent() { this->Save(true); }
	};

    class OutputGeneratorException : public DREAM::FVM::FVMException {
    public:
        template<typename ... Args>
        OutputGeneratorException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("OutputGenerator");
        }
    };
}

#endif/*_DREAM_OUTPUT_GENERATOR_HPP*/
