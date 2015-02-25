#ifndef _MODELFITTER_H_
#define _MODELFITTER_H_

#include <string>
#include <vector>

#include "libsrc/MathMatrix.h"
#include "Result.h"
#include "Summary.h"

class DataConsolidator;
// take X, Y, Cov and fit model
// note, ModelFitter will use VCFData as READ-ONLY data structure,
// and collapsing results are stored internally.
class ModelFitter{
 public:
  virtual int fit(DataConsolidator* dc) = 0;

  // write result header
  virtual void writeHeader(FileWriter* fp, const Result& siteInfo) = 0;
  // write model output
  virtual void writeOutput(FileWriter* fp, const Result& siteInfo) = 0;
  // write footnotes in model output
  virtual void writeFootnote(FileWriter* fp) {};

  ModelFitter(){
    this->modelName = "UninitializedModel";
    this->binaryOutcome = false; // default: using continuous outcome
    this->familyModel = false;   // is this case for for family only?
    this->indexResult = false;
  }
  virtual ~ModelFitter() {}
  const std::string& getModelName() const { return this->modelName; };
  // for particular class to call when fitting repeatedly
  // e.g. clear permutation counter
  // e.g. clear internal cache
  virtual void reset() { this->result.clearValue();};
  // virtual void needFittingCovariate();
  bool isBinaryOutcome() const{
    return this->binaryOutcome;
  }
  void setBinaryOutcome() {
    this->binaryOutcome = true;
  }
  bool isFamilyModel() const {
    return this->familyModel;
  }
  bool needToIndexResult() const {
    return this->indexResult;
  }
  void setContinuousOutcome() {
    this->binaryOutcome = false;
  }
  const Result& getResult() const {
    return this->result;
  }
  void appendHeader(SummaryHeader* h) {
    header.push_back(h);
  }
  void setPrefix(const std::string& p) {
    this->outputPrefix = p;
  }
 protected:
  std::string modelName;
  bool binaryOutcome;
  bool familyModel;
  bool indexResult;
  Result result;
  FileWriter* fp;

  // optionally used
  std::string outputPrefix;
  std::vector<SummaryHeader*> header;
}; // end ModelFitter
#endif /* _MODELFITTER_H_ */
