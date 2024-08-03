#ifndef _DREAM_OUTPUT_GENERATOR_SFILE_HPP
#define _DREAM_OUTPUT_GENERATOR_SFILE_HPP

#include <softlib/SFile.h>
#include <string>
#include "DREAM/OutputGenerator.hpp"

namespace DREAM {
	class OutputGeneratorSFile : public OutputGenerator {
	protected:
        std::string filename;
		SFile *sf=nullptr;

		virtual void SaveGrids(const std::string&, bool) override;
		virtual void SaveIonMetaData(const std::string&) override;
		virtual void SaveOtherQuantities(const std::string&) override;
		virtual void SaveSettings(const std::string&) override;
        virtual void SaveSolverData(const std::string&) override;
		virtual void SaveTimings(const std::string&) override;
		virtual void SaveUnknowns(const std::string&, bool) override;

		void SaveEquilibrium(SFile*, const std::string&);
        void SaveMomentumGrid(SFile*, const std::string&, FVM::Grid*, enum OptionConstants::momentumgrid_type);
        void WriteCopyArray(SFile*, const std::string&, const real_t *const*, const len_t, const len_t);
        void WriteCopyMultiArray(SFile*, const std::string&, const real_t *const*, const sfilesize_t, const sfilesize_t[]);
	public:
		OutputGeneratorSFile(EquationSystem*, const std::string&, bool savesettings=true);
		OutputGeneratorSFile(EquationSystem*, SFile*, bool savesettings=true);
        virtual ~OutputGeneratorSFile();

        SFile *GetSFile() { return this->sf; }

        virtual void Save(bool current=false) override;
	};
}

#endif/*_DREAM_OUTPUT_GENERATOR_SFILE_HPP*/
