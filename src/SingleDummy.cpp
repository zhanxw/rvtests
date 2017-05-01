#include "SingleDummy.h"

#include "src/DataConsolidator.h"

// this function fits the model
int SingleDummy::fit(DataConsolidator* dc) {
  // this is a simple counter
  counter_++;

  // this is extract genotype, phenotype and covariate for further analysis
  // let's say you want to get the max and min values among them
  Matrix& phenotype = dc->getPhenotype();
  Matrix& genotype = dc->getGenotype();
  Matrix& cov = dc->getCovariate();
  max_ = std::max(std::max(phenotype.Max(), genotype.Max()), cov.Max());
  min_ = std::min(std::max(phenotype.Min(), genotype.Min()), cov.Min());

  return 0;
}

// write result header
void SingleDummy::writeHeader(FileWriter* fp, const Result& siteInfo) {
  siteInfo.writeHeaderTab(fp);
  result.writeHeaderLine(fp);
}
// write model output
void SingleDummy::writeOutput(FileWriter* fp, const Result& siteInfo) {
  // write site info
  siteInfo.writeValueTab(fp);

  // write model results
  result.updateValue("NCOUNT", counter_);
  result.updateValue("MAX", max_);
  result.updateValue("MIN", min_);
  result.writeValueLine(fp);
}

// optional functions //
// write footnotes in model output
// virtual void writeFootnote(FileWriter* fp);

//   // if you need to parse parameters, change here
// virtual int setParameter(const ModelParser& parser) {
//     // maybe check whether user provides duplicated/non-existing parameters
//     return 0;
//   }

SingleDummy::SingleDummy() {
  counter_ = 0;
  result.addHeader("NCOUNT");
  result.addHeader("MAX");
  result.addHeader("MIN");
}
