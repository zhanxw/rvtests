#ifndef _MODELFITTER_H_
#define _MODELFITTER_H_

/* #include <gsl/gsl_rng.h> */
/* #include <gsl/gsl_randist.h> */

/* #include <deque> */

#include <string>
#include <vector>

#include "libsrc/MathMatrix.h"
#include "Result.h"
#include "Summary.h"

#if 0
#include "base/ParRegion.h"
#include "regression/LogisticRegression.h"
#include "regression/LogisticRegressionScoreTest.h"
#include "regression/LogisticRegressionVT.h"
#include "regression/LinearRegression.h"
#include "regression/LinearRegressionScoreTest.h"
#include "regression/LinearRegressionVT.h"
#include "regression/Skat.h"
#include "regression/Table2by2.h"
#include "regression/kbac_interface.h"
#include "regression/FastLMM.h"
#include "regression/GrammarGamma.h"
#include "regression/MetaCov.h"
#include "regression/MatrixOperation.h"
#include "regression/MultivariateVT.h"
#include "regression/FamSkat.h"
#include "regression/FirthRegression.h"

#include "DataConsolidator.h"
#include "LinearAlgebra.h"
#include "ModelUtil.h"
#include "snp_hwe.c"
#include "Result.h"
#include "Summary.h"
#endif

#if 0
// may decrease speed.
#ifdef _OPENMP
#include <omp.h>
#pragma message "Enable multithread using OpenMP"
#endif
#endif

#if 0
extern SummaryHeader* g_SummaryHeader;

typedef void (*CollapsingFunction)(Matrix& in, const std::vector<int>& idx, Matrix* out, int index);

// various collapsing method
// they all take people by marker matrix
// and they won't take special care of missing genotypes
double getMarkerFrequency(Matrix& in, int col);
void getMarkerFrequency(Matrix& in, std::vector<double>* freq);
double getMarkerFrequencyFromControl(Matrix& in, Vector& pheno, int col);

void cmcCollapse(Matrix& in, Matrix* out);
void cmcCollapse(Matrix& in, const std::vector<int>& idx, Matrix* out, int index);

void zegginiCollapse(Matrix& in, Matrix* out);
void zegginiCollapse(Matrix& in, const std::vector<int>& idx, Matrix* out, int index);

void fpCollapse(Matrix& in, Matrix* out);

void madsonBrowningCollapse(Matrix& genotype, Vector& phenotype, Matrix* out);

void groupFrequency(const std::vector<double>& freq, std::map<double, std::vector<int> >* group);
void convertToMinorAlleleCount(Matrix& in, Matrix* g);
void convertToReferenceAlleleCount(Matrix& in, Matrix* g);

void makeVariableThreshodlGenotype(Matrix& in,
                                   const std::vector<double>& freqIn,
                                   Matrix* out,
                                   std::vector<double>* freqOut,
                                   void (*collapseFunc)(Matrix& , const std::vector<int>& , Matrix*, int)
                                   );
void appendHeritability(FileWriter* fp, const FastLMM& model);
void appendHeritability(FileWriter* fp, const GrammarGamma& model);
#endif

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
