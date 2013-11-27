#include "VCFInputFile.h"
#include "VCFFilter.h"
#include <string>

class VCFExtractor: public VCFInputFile, public VCFSiteFilter {
 public:
  VCFExtractor(const std::string& fn): VCFInputFile(fn){
  };
  virtual ~VCFExtractor() {};
  bool passFilter();
};
