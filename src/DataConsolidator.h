#ifndef _DATACONSOLIDATOR_H_
#define _DATACONSOLIDATOR_H_

#include "base/KinshipHolder.h"
#include "base/Logger.h"
#include "base/MathMatrix.h"
#include "base/ParRegion.h"
#include "libsrc/Random.h"
#include "regression/Formula.h"
#include "src/GenotypeCounter.h"
#include "src/Result.h"

class SimpleMatrix;

extern Logger* logger;

class EigenMatrix;

/**
 * @return true if any of the markers (@param col) of @param genotype (people by
 * marker) is missing
 */
inline bool hasMissingMarker(Matrix& genotype, int col);

/**
 * Remove columns of markers in @param genotype (people by marker) where there
 * are missing genotypes
 */
void removeMissingMarker(Matrix* genotype);

/**
 * @return true if markers on @param col of @param genotype (people by marker)
 * is monomorphic (genotypes are all the same)
 */
bool isMonomorphicMarker(Matrix& genotype, int col);

/**
 * remove monomorphic columns of @param genotype
 */
void removeMonomorphicMarker(Matrix* genotype);

/**
 * Convert genotype back to reference allele count
 * e.g. genotype 2 means homAlt/homAlt, so it has reference allele count 0
 */
void convertToMinorAlleleCount(Matrix& in, Matrix* g);

/**
 * This class is in charge of cleanning data before fitting in model
 * The cleaning step includes:
 *  remove monomorphic sites
 *  handle missing data genotype (impute to mean, impute by HWE, or
 *   filter out mssing genotypes and its corresponding phenotypes, covariates)
 *  (future todo) include weights (GERP, Sift)
 *  (future todo) re-weight genotype (dominate model, recessive model)
 */
class DataConsolidator {
 public:
  const static int UNINITIALIZED = 0;
  const static int IMPUTE_MEAN = 1;
  const static int IMPUTE_HWE = 2;
  const static int DROP = 3;
  const static int KINSHIP_AUTO = 0;
  const static int KINSHIP_X = 1;
  typedef enum { ANY_SEX = -1, MALE = 1, FEMALE = 2 } PLINK_SEX;
  typedef enum { ANY_PHENO = -1, CTRL = 1, CASE = 2 } PLINK_PHENOTYPE;

 public:
  DataConsolidator();
  virtual ~DataConsolidator();
  void setStrategy(const int s) { this->strategy = s; };

  /**
   * Impute missing genotype (<0) according to population frequency (p^2, 2pq,
   * q^2)
   */
  void imputeGenotypeByFrequency(Matrix* genotype, Random* r);

  /**
   * Impute missing genotype (<0) according to its mean genotype
   * @param genotype (people by marker matrix)
   */
  void imputeGenotypeToMean(Matrix* genotype);

  /**
   * @param pheno, @param cov @param genotype are all ordered and sorted by the
   * same set of samples
   *
   * NOTE: we assume @param pheno, @param cov are not changed outside of this
   * function;
   * we assume @param geno is always changed
   */
  void consolidate(Matrix& pheno, Matrix& cov, Matrix& geno);
  /**
      * Compare @param a and @param b by comparing their common finite elements.
      */
  bool isEqual(Matrix& a, Matrix& b);
  /**
   * @param g , row @param r: all elements are non-missing
   */
  bool isNoMissingGenotypeInRow(Matrix& g, int r) {
    const int n = g.cols;
    for (int i = 0; i < n; ++i) {
      if (g(r, i) < 0) return false;
    }
    return true;
  }
  /**
   * Copy the row @param srcRow in matrix @param src to @param dst at row
   * @destRow
   * NOTE: @param dest rows and columns may be changed depends on the dimension
   * of @param src and @param srcRow
   */
  void copyRow(Matrix& src, const int srcRow, Matrix* dest, const int destRow);

  void copyColName(Matrix& src, Matrix* dest) {
    dest->Dimension(dest->rows, src.cols);
    for (int i = 0; i < src.cols; ++i) {
      dest->SetColumnLabel(i, src.GetColumnLabel(i));
    }
  }
  void setPhenotypeName(const std::vector<std::string>& name) {
    this->originalRowLabel = name;
    this->rowLabel = name;
  }
  const std::vector<std::string>& getRowLabel() const { return this->rowLabel; }
  Matrix& getGenotype() { return this->genotype; }
  Matrix& getFlippedToMinorPolymorphicGenotype() {
    convertToMinorAlleleCount(this->genotype, &this->flippedToMinorGenotype);
    removeMonomorphicMarker(&flippedToMinorGenotype);
    return this->flippedToMinorGenotype;
  }
  Matrix& getOriginalGenotype() { return this->originalGenotype; }
  Matrix& getPhenotype() { return this->phenotype; }
  Matrix& getCovariate() { return this->covariate; }
  Vector& getWeight() { return this->weight; }
  Result& getResult() { return this->result; }

  /**
   * Count @param homRef, @param het, @param homAlt and @param missing
   * from the genotype column specified by @param columnIndex
   * @param sex : only process specified sex (1, male; 2, female; <0, any)
   * @param phenotype: only process specified phenotype (1, control; 2, case; <0
   * any)
   * @return 0 if succeed
   */
  int countRawGenotype(int columnIndex, const PLINK_SEX sex,
                       const PLINK_PHENOTYPE phenotype,
                       GenotypeCounter* counter) const;

  int countRawGenotypeFromCase(int columnIndex,
                               GenotypeCounter* counter) const {
    return countRawGenotype(columnIndex,
                            ANY_SEX,  // PLINK male
                            CASE, counter);
  }

  int countRawGenotypeFromControl(int columnIndex,
                                  GenotypeCounter* counter) const {
    return countRawGenotype(columnIndex, ANY_SEX,
                            CTRL,  // any phenotype
                            counter);
  }

  int countRawGenotype(int columnIndex, GenotypeCounter* counter) const {
    return countRawGenotype(columnIndex, ANY_SEX,
                            ANY_PHENO,  // any phenotype
                            counter);
  }
  int countRawGenotypeFromFemale(int columnIndex,
                                 GenotypeCounter* counter) const {
    return countRawGenotype(columnIndex, FEMALE,
                            ANY_PHENO,  // any phenotype
                            counter);
  }
  int countRawGenotypeFromFemaleCase(int columnIndex,
                                     GenotypeCounter* counter) const {
    return countRawGenotype(columnIndex, FEMALE, CASE, counter);
  }
  int countRawGenotypeFromFemaleControl(int columnIndex,
                                        GenotypeCounter* counter) const {
    return countRawGenotype(columnIndex, FEMALE, CTRL, counter);
  }

  bool isPhenotypeUpdated() const { return this->phenotypeUpdated; }
  bool isCovariateUpdated() const { return this->covariateUpdated; }

  void setParRegion(ParRegion* p) { this->parRegion = p; }

  // Sex (1=male; 2=female; other=unknown)
  void setSex(const std::vector<int>* sex) { this->sex = sex; };
  void setFormula(const FormulaVector* formula) { this->formula = formula; };
  const FormulaVector* getFormula() const { return this->formula; };
  void setGenotypeCounter(const std::vector<GenotypeCounter>& c) {
    this->counter = &c;
  };
  /**
   * Check if genotype matrix column @param columnIndex is a chromosome X.
   */
  bool isHemiRegion(int columnIndex) {
    assert(this->parRegion);
    std::string chromPos = this->genotype.GetColumnLabel(columnIndex);
    size_t posColon = chromPos.find(":");
    if (posColon == std::string::npos) return false;
    std::string chrom = chromPos.substr(0, posColon);
    int pos = atoi(chromPos.substr(posColon + 1));
    return this->parRegion->isHemiRegion(chrom, pos);
  }
  /**
   * Recode this->orginalGenotype to @parma geno, using pre-specified imputation
   * strategy
   */
  void codeGenotypeForDominantModel(Matrix* geno);
  void codeGenotypeForRecessiveModel(Matrix* geno);

 public:
  // codes to check before regression
  int preRegressionCheck(Matrix& pheno, Matrix& cov);
  int checkColinearity(Matrix& cov);
  int checkPredictor(Matrix& pheno, Matrix& cov);

 public:
  double getMarkerFrequency(int col);
  void getMarkerFrequency(std::vector<double>* freq);

 public:
  // codes related to kinship
  int setKinshipSample(const std::vector<std::string>& samples);
  int setKinshipFile(int kinshipType, const std::string& fileName);
  int setKinshipEigenFile(int kinshipType, const std::string& fileName);
  int loadKinship(int kinshipType);

  const EigenMatrix* getKinshipForAuto() const {
    return this->kinship[KINSHIP_AUTO].getK();
  }
  const EigenMatrix* getKinshipUForAuto() const {
    return this->kinship[KINSHIP_AUTO].getU();
  }
  const EigenMatrix* getKinshipSForAuto() const {
    return this->kinship[KINSHIP_AUTO].getS();
  }
  bool hasKinshipForAuto() const {
    return this->kinship[KINSHIP_AUTO].isLoaded();
  }

  const EigenMatrix* getKinshipForX() const {
    return this->kinship[KINSHIP_X].getK();
  }
  const EigenMatrix* getKinshipUForX() const {
    return this->kinship[KINSHIP_X].getU();
  }
  const EigenMatrix* getKinshipSForX() const {
    return this->kinship[KINSHIP_X].getS();
  }
  bool hasKinshipForX() const { return this->kinship[KINSHIP_X].isLoaded(); }

  bool hasKinship() const {
    return this->hasKinshipForAuto() || this->hasKinshipForX();
  }

 public:
  /**
   * Create data files for BOLT-LMM heritability estimation, inlcudes
   *  - a set of binary PLINK file (genotype, phentoype)
   *  - an optional .covar file that stores covariates
   * the outputted sample has to follow the given order in @param sampleName
   */
  int prepareBoltModel(const std::string& prefix,
                       const std::vector<std::string>& sampleName,
                       const SimpleMatrix& phenotype);
  const std::string& getBoltGenotypeFilePrefix() const {
    return this->boltPrefix;
  }
#if 0
  /**
   * Load sample by genotype matrix
   */
  int loadGenotype(const std::string& prefix);
  /**
   * Load sample by genotype matrix, fill missing to mean, and equalize variance
   */
  int loadNormalizedGenotype(const std::string& prefix);
  EigenMatrix* getFullGenotype();
#endif

 private:
  // don't copy
  DataConsolidator(const DataConsolidator&);
  DataConsolidator& operator=(const DataConsolidator&);

 private:
  int strategy;
  Random random;
  Matrix genotype;
  Matrix flippedToMinorGenotype;
  Matrix phenotype;
  Matrix covariate;
  Vector weight;
  Result result;
  Matrix originalGenotype;
  bool phenotypeUpdated;
  bool covariateUpdated;
  std::vector<std::string> originalRowLabel;
  std::vector<std::string> rowLabel;
  KinshipHolder kinship[2];  // 2: include both AUTO and X kinships
  std::string boltPrefix;    // prefix for a set of file for BoltLMM model

  // sex chromosome adjustment
  const std::vector<int>* sex;
  // store formulae
  const FormulaVector* formula;
  const std::vector<GenotypeCounter>* counter;
  ParRegion* parRegion;
};  // end DataConsolidator

#endif /* _DATACONSOLIDATOR_H_ */
