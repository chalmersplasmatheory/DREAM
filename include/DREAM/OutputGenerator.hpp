#ifndef _DREAM_OUTPUT_GENERATOR_HPP
#define _DREAM_OUTPUT_GENERATOR_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/OtherQuantityHandler.hpp"

namespace DREAM {
	class OutputGenerator {
	protected:
		FVM::Grid *grid;
		FVM::UnknownQuantityHandler *unknowns;
		IonHandler *ions;
		OtherQuantityHandler *oqty;
		EquationSystem *eqsys;

		virtual void SaveGrids(const std::string&) = 0;
		virtual void SaveIonMetaData(const std::string&) = 0;
		virtual void SaveOtherQuantities(const std::string&) = 0;
		virtual void SaveSettings(const std::string&) = 0;
		virtual void SaveTimings(const std::string&) = 0;
		virtual void SaveUnknowns(const std::string&) = 0;
	public:
		OutputGenerator(
			FVM::Grid*, FVM::UnknownQuantityHandler*,
			IonHandler*, OtherQuantityHandler*,
			EquationSystem*
		);

		virtual void Save();
	};
}

#endif/*_DREAM_OUTPUT_GENERATOR_HPP*/
