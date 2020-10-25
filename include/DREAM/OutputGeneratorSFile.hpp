#ifndef _DREAM_OUTPUT_GENERATOR_SFILE_HPP
#define _DREAM_OUTPUT_GENERATOR_SFILE_HPP

#include <softlib/SFile.h>
#include "DREAM/OutputGenerator.hpp"

namespace DREAM {
	class OutputGeneratorSFile : public OutputGenerator {
	protected:
		SFile *sf;

		virtual void SaveGrids(const std::string&);
		virtual void SaveIonMetaData(const std::string&);
		virtual void SaveOtherQuantities(const std::string&);
		virtual void SaveSettings(const std::string&);
		virtual void SaveTimings(const std::string&);
		virtual void SaveUnknowns(const std::string&);
	public:
		OutputGeneratorSFile(
			FVM::Grid*, FVM::UnknownQuantityHandler*,
			IonHandler*, OtherQuantityHandler*,
			EquationSystem*
		);
	};
}

#endif/*_DREAM_OUTPUT_GENERATOR_SFILE_HPP*/
