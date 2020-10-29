#ifndef _DREAM_OUTPUT_GENERATOR_SFILE_HPP
#define _DREAM_OUTPUT_GENERATOR_SFILE_HPP

#include <softlib/SFile.h>
#include "DREAM/OutputGenerator.hpp"

namespace DREAM {
	class OutputGeneratorSFile : public OutputGenerator {
	protected:
		SFile *sf;

		virtual void SaveGrids(const std::string&, bool);
		virtual void SaveIonMetaData(const std::string&);
		virtual void SaveOtherQuantities(const std::string&);
		virtual void SaveSettings(const std::string&);
		virtual void SaveTimings(const std::string&);
		virtual void SaveUnknowns(const std::string&, bool);

        void SaveMomentumGrid(SFile*, const std::string&, FVM::Grid*, enum OptionConstants::momentumgrid_type);
        void WriteCopyArray(SFile*, const std::string&, const real_t *const*, const len_t, const len_t);
        void WriteCopyMultiArray(SFile*, const std::string&, const real_t *const*, const sfilesize_t, const sfilesize_t[]);
	public:
		OutputGeneratorSFile(EquationSystem*, const std::string&);
		OutputGeneratorSFile(EquationSystem*, SFile*);
        virtual ~OutputGeneratorSFile();
	};
}

#endif/*_DREAM_OUTPUT_GENERATOR_SFILE_HPP*/
