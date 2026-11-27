#ifndef _DREAM_STATIC_QUANTITY_HPP
#define _DREAM_STATIC_QUANTITY_HPP

#include <functional>
#include <string>
#include <softlib/SFile.h>


namespace DREAM {
	class StaticQuantity {
	protected:
		std::string name, description;
		std::function<void(SFile*, const std::string&)> writeFunction;

		bool active = false;

	public:
		StaticQuantity(
			const std::string& name, 
			const std::string& description, 
			std::function<void(SFile*, const std::string&)> writeFunction
		) : name(name), description(description) {
			this->writeFunction = writeFunction;
		}

		~StaticQuantity() {}

		void Activate() { this->active = true; }
		bool IsActive() { return this->active; }
		const std::string& GetName() { return this->name; }
		const std::string& GetDescription() { return this->description; }

		void WriteSFile(SFile *sf, const std::string& path) {
			this->writeFunction(sf, path +  this->name);
		}
	};
}

#endif/*_DREAM_STATIC_QUANTITY_HPP*/
