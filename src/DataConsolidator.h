#ifndef _DATACONSOLIDATOR_H_
#define _DATACONSOLIDATOR_H_

#include "base/Logger.h"
#include "base/ParRegion.h"
#include "libsrc/MathMatrix.h"
#include "libsrc/Random.h"

#include "base/KinshipHolder.h"
#include "regression/Formula.h"
#include "src/GenotypeCounter.h"
#include "src/Result.h"

class SimpleMatrix;

extern Logger* logger;

class WarningOnce {
 public:
  WarningOnce(const std::string& msg) : warningGiven(false), msg(msg){};
  void warningIf(bool cond) {
    if (cond && !warningGiven) {
      warningGiven = true;
      fprintf(stderr, "%s", msg.c_str());
    }
  }

 private:
  bool warningGiven;
  std::string msg;
};

class EigenMatrix;

/**
 * @return true if any of the markers (@param col) of @param genotype (people by
 * marker) is missing
 */
inline bool hasMissingMarker(Matrix& genotype, int col) {
  if (col >= genotype.cols || col < 0) {
    logger->error("Invalid check of missing marker.");
    return false;
  }

  for (int r = 0; r < genotype.rows; ++r) {
    if (genotype[r][col] < 0) return true;
  }
  return false;
};

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
  void consolidate(Matrix& pheno, Matrix& cov, Matrix& geno) {
    if (&geno != &this->originalGenotype) {
      this->originalGenotype = geno;
      copyColName(geno, &this->originalGenotype);
      // fprintf(stderr, "== Copy occured\n");
    }

    this->genotype = geno;
    copyColName(geno, &this->genotype);

    if (isPhenotypeUpdated()) {
      copyColName(pheno, &this->phenotype);
    }
    if (isCovariateUpdated()) {
      copyColName(cov, &this->covariate);
    }

    // impute missing genotypes
    if (this->strategy == IMPUTE_MEAN) {
      // impute missing genotypes
      imputeGenotypeToMean(&this->genotype);

      // handle phenotype
      if (isPhenotypeUpdated()) {
        this->phenotypeUpdated = !isEqual(this->phenotype, pheno);
        this->phenotype = pheno;
      } else {
        // no need to update phenotype
      }

      // handle covariate
      if (isCovariateUpdated()) {
        this->covariateUpdated = !isEqual(this->covariate, cov);
        this->covariate = cov;
      } else {
        // no need to update covariate
      }
    } else if (this->strategy == IMPUTE_HWE) {
      // impute missing genotypes
      imputeGenotypeByFrequency(&genotype, &this->random);
      // handle phenotype
      if (isPhenotypeUpdated()) {
        this->phenotypeUpdated = !isEqual(this->phenotype, pheno);
        this->phenotype = pheno;
      } else {
        // no need to update phenotype
      }

      // handle covariate
      if (isCovariateUpdated()) {
        this->covariateUpdated = !isEqual(this->covariate, cov);
        this->covariate = cov;
      } else {
        // no need to update covariate
      }
    } else if (this->strategy == DROP) {
      // (TODO) should also consider how kinship matrix changes.

      // we process genotype matrix (people by marker)
      // if for the same people, any marker is empty, we will remove this people
      int idxToCopy = 0;
      for (int i = 0; i < (genotype).rows; ++i) {
        if (isNoMissingGenotypeInRow(genotype, i)) {
          copyRow(genotype, i, &genotype, idxToCopy);
          copyRow(cov, i, &covariate, idxToCopy);
          copyRow(pheno, i, &phenotype, idxToCopy);
          rowLabel[idxToCopy] = originalRowLabel[i];
          idxToCopy++;
        }
      }
      genotype.Dimension(idxToCopy, genotype.cols);
      covariate.Dimension(idxToCopy, cov.cols);
      phenotype.Dimension(idxToCopy, pheno.cols);
      this->phenotypeUpdated = true;
      this->covariateUpdated = true;
    } else {
      logger->error(
          "Uninitialized consolidation methods to handle missing data!");
    }
  }  // end consolidate
     /**
      * Compare @param a and @param b by comparing their common finite elements.
      */
  bool isEqual(Matrix& a, Matrix& b) {
    if (a.rows != b.rows) return false;
    if (a.cols != b.cols) return false;
    const int nr = a.rows;
    const int nc = a.cols;
    for (int i = 0; i < nr; ++i) {
      for (int j = 0; j < nc; ++j) {
        if (finite(a[i][j]) && finite(b[i][j]) && a[i][j] != b[i][j])
          return false;
      }
    }
    return true;
  }
  /**
   * @param g , row @param r: all elements are non-missing
   */
  bool isNoMissingGenotypeInRow(Matrix& g, int r) {
    const int n = g.cols;
    for (int i = 0; i < n; ++i) {
      if (g[r][i] < 0) return false;
    }
    return true;
  }
  void copyRow(Matrix& src, const int srcRow, Matrix* dest, const int destRow) {
    Matrix& m = *dest;
    if (m.cols < src.cols) {
      m.Dimension(m.rows, src.cols);
    }
    if (m.rows <= destRow) {
      m.Dimension(destRow + 1, m.cols);
    }
    const int n = m.cols;
    for (int i = 0; i < n; ++i) {
      m[destRow][i] = src[srcRow][i];
    }
  }
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
                       GenotypeCounter* counter) const {
    if (columnIndex < 0 || columnIndex >= originalGenotype.cols) {
      return -1;
    }
    if (sex > 0 && sex != MALE && sex != FEMALE) return -2;
    if (sex > 0 && (int)this->sex->size() != originalGenotype.rows) return -3;
    if (phenotype > 0 && phenotype != CTRL && phenotype != CASE) return -2;

    for (int i = 0; i < originalGenotype.rows; ++i) {
      if (sex > 0 && (*this->sex)[i] != sex) {
        continue;
      }
      // + 1: PLINK use 1 and 2 as ctrl and case, but
      // internally, we use 0 and 1.
      if (phenotype > 0 && (int)(this->phenotype[i][0] + 1) != phenotype) {
        continue;
      }

      const double& g = originalGenotype[i][columnIndex];
      counter->add(g);
    }

    return 0;  // success
  }

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

  void codeGenotypeForDominantModel(Matrix* geno) {
    int n = genotype.cols;
    static WarningOnce warning("Encoding only use the first variant!\n");
    warning.warningIf(n != 1);

    int m = genotype.rows;
    if (n != 1) {
      fprintf(stderr, "n = %d, m = %d \n", n, m);
    }

    geno->Dimension(m, 1);
    double s = 0;  // sum of genotypes
    int numGeno = 0;
    if (this->strategy == IMPUTE_MEAN || this->strategy == IMPUTE_HWE) {
      for (int i = 0; i < m; ++i) {
        if (this->originalGenotype[i][0] < 0) continue;

        if (this->originalGenotype[i][0] > 0.5) {
          (*geno)[i][0] = 1.0;
          s += 1.;
          numGeno++;
        } else {
          (*geno)[i][0] = 0.0;
          numGeno++;
        }
      }
      double avg = 0.0;
      if (numGeno > 0) {
        avg = s / numGeno;
      }
      for (int i = 0; i < m; ++i) {
        if (this->originalGenotype[i][0] < 0) (*geno)[i][0] = avg;
      }
    } else if (this->strategy == DROP) {
      for (int i = 0; i < m; ++i) {
        if (this->genotype[i][0] > 0.5) {
          (*geno)[0][i] = 1.;
        } else {
          (*geno)[0][i] = 0.;
        }
      }
    }
  }
  void codeGenotypeForRecessiveModel(Matrix* geno) {
    int n = genotype.cols;
    static WarningOnce warning("Encoding only use the first variant!\n");
    warning.warningIf(n != 1);

    int m = genotype.rows;
    geno->Dimension(m, 1);
    double s = 0;  // sum of genotypes
    int numGeno = 0;
    if (this->strategy == IMPUTE_MEAN || this->strategy == IMPUTE_HWE) {
      for (int i = 0; i < m; ++i) {
        if (this->originalGenotype[i][0] < 0) continue;

        if (this->originalGenotype[i][0] > 1.5) {
          (*geno)[i][0] = 1.0;
          s += 1.;
          numGeno++;
        } else {
          (*geno)[i][0] = 0.0;
          numGeno++;
        }
      }
      double avg = 0.0;
      if (numGeno > 0) {
        avg = s / numGeno;
      }
      for (int i = 0; i < m; ++i) {
        if (this->originalGenotype[i][0] < 0) (*geno)[i][0] = avg;
      }
    } else if (this->strategy == DROP) {
      for (int i = 0; i < m; ++i) {
        if (this->genotype[i][0] > 1.5) {
          (*geno)[0][i] = 1.;
        } else {
          (*geno)[0][i] = 0.;
        }
      }
    }
  }

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
