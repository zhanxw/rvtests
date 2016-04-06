#ifndef _DATALOADER_H_
#define _DATALOADER_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#include "regression/Formula.h"
#include "base/SimpleMatrix.h"
#include "libsrc/MathVector.h"
#include "libsrc/MathMatrix.h"

class DataLoader {
 public:
  typedef enum {
    COVARIATE_IMPUTE,
    COVARIATE_DROP,
    COVARIATE_KEEP
  } HandleMissingCov;
  typedef enum {
    PHENOTYPE_QTL,
    PHENOTYPE_BINARY,
    PHENOTYPE_UNKNOWN
  } PhenotypeType;

 public:
  DataLoader();
  // phenotypes related
  int loadPhenotype(const std::string& pheno, const std::string& mpheno,
                    const std::string& phenoName);
  int arrangePhenotype(const std::vector<std::string>& names,
                       std::vector<std::string>* droppedNames);

  // covariates related
  void setImputeCovariate();
  int loadCovariate(const std::string& covar, const std::string& covName);
  int arrangeCovariate(const std::vector<std::string>& names,
                       std::vector<std::string>* droppedNames);

  int loadSex();
  int useSexAsCovariate();
  int loadMarkerAsCovariate(const std::string& inVcf,
                            const std::string& marker);

  // load multiple phenotype
  int loadMultiplePhenotype(const std::string& multiplePhenotype,
                            const std::string& pheno, const std::string& covar);

  // sanity check
  int checkConstantCovariate();

  // transformations
  int useResidualAsPhenotype();
  int inverseNormalizePhenotype();

  // phenotype-related utilities
  PhenotypeType detectPhenotypeType() const;
  int setTraitType(PhenotypeType t);
  bool isBinaryPhenotype() const { return binaryPhenotype; };

  // setters
  void setPhenotypeImputation(bool b);
  void setCovariateImputation(bool b);

  // getters
  const SimpleMatrix& getPhenotype() { return this->phenotype; };
  const SimpleMatrix& getCovariate() { return this->covariate; };
  const std::vector<int>& getSex() { return this->sex; };
  const FormulaVector& getFormula() const { return this->formula; };

 private:
  // for multiple traits
  std::vector<SimpleMatrix> phenotypes;  // sample by traits
  std::vector<SimpleMatrix> covariates;  // sample by covariates
  FormulaVector formula;

  SimpleMatrix& phenotype;  // sample by traits
  SimpleMatrix& covariate;  // sample by covariates
  bool binaryPhenotype;
  std::vector<int> sex;  // plink coded genders

  // external parameters
  std::string FLAG_pheno;
  std::string FLAG_mpheno;
  std::string FLAG_phenoName;
  bool FLAG_imputePheno;
  std::string FLAG_multiplePheno;
  
  std::string FLAG_cov;
  std::string FLAG_covName;
  bool FLAG_imputeCov;

  std::string FLAG_inVcf;
  std::string FLAG_condition;

  // intermediate values
  std::set<std::string> sampleToDropInCovariate;
};

/**
 * Extract covaraite from file @param fn.
 * Only samples included in @param includedSample will be processed
 * If some samples appear more than once, only the first appearance will be
 * readed
 * Only covaraites provided in @param covNameToUse will be included
 * Missing values will be imputed to the mean columnwise.
 * Result will be put to @param mat (sample by covariate) and @param
 * sampleToDrop
 * @return number of sample loaded (>=0); or a minus number meaning error
 * @param sampleToDrop: store samples that are not found in covariate.
 */
int extractCovariate(const std::string& fn,
                     const std::vector<std::string>& sampleToInclude,
                     const std::vector<std::string>& covNameToUse,
                     DataLoader::HandleMissingCov handleMissingCov,
                     SimpleMatrix* mat, std::set<std::string>* sampleToDrop);

/**
 * Load covariate from @param fn, using specified @param covNameToUse, for given
 * @param includedSample
 * covariate will be stored in @param covariate, and column names will be stored
 * in @colNames
 * if covariate file missed some samples, those sample names will be stored in
 * @sampleToDrop
 * NOTE: for missing values in a covariate, it will drop this covariate out of
 * the following anaylysis
 * @return number of samples have covariates.
 * Example:
 * includedSample = [A, B, C] and in covariate file we have [B, C, C, D]
 * then output covariate have 3 rows corresponding to [A, B, C]
 * row C filled by the last C in covariate file
 * sample D will be in sampleToDrop
 */
int _loadCovariate(const std::string& fn,
                   const std::vector<std::string>& includedSample,
                   const std::vector<std::string>& covNameToUse,
                   DataLoader::HandleMissingCov handleMissingCov,
                   SimpleMatrix* covariate, std::vector<std::string>* colNames,
                   std::set<std::string>* sampleToDrop);

int _loadCovariate(const std::string& fn,
                   const std::vector<std::string>& includedSample,
                   const std::string& covNameToUse,
                   DataLoader::HandleMissingCov handleMissingCov,
                   SimpleMatrix* covariate, std::vector<std::string>* colNames,
                   std::set<std::string>* sampleToDrop);

/**
 * @return number of phenotypes read. -1 if errors
 * @param phenoCol, which phenotype column to use, similar to plink, it should
 * be the order of phenotype, e.g. phenoCol = 2, meaning the second phenotype
 * @param phenoName, which phenotype header to use.
 */
int loadPedPhenotypeByColumn(const char* fn, std::map<std::string, double>* p,
                             int phenoCol);

/**
 * @return number of phenotypes read. -1 if errors
 * @param phenoName, which phenotype header to use.
 *
 */
int loadPedPhenotypeByHeader(const char* fn, std::map<std::string, double>* p,
                             const char* phenoHeader);

/**
 * @return true if @param phenotype is either:  1: unaffected, 2: affected,  -9,
 * 0: missing
 */
bool _isBinaryPhenotype(const std::vector<double>& phenotype);

/**
 * Convert binary phenotype 1,2 (PLINK format) to 0,1 (logistic regression)
 */
bool convertBinaryPhenotype(std::vector<double>* p);

/**
 * according to the order of @param vcfSampleNames, put phenotypes to @param
 * phenotypeInOrder
 * @param imputePhenotype: if true, we will impute phenotpye to the average for
 * those have genotype but no phenotype;
 *                         if false, we will drop those samples
 */
void rearrange(const std::map<std::string, double>& phenotype,
               const std::vector<std::string>& vcfSampleNames,
               std::vector<std::string>* vcfSampleToDrop,
               std::vector<std::string>* phenotypeNameInOrder,
               std::vector<double>* phenotypeValueInOrder,
               bool imputePhenotype);

int _loadSex(const std::string& fn,
             const std::vector<std::string>& includedSample,
             std::vector<int>* sex);

/**
 * when @param sex does not equal to 1 (male) or 2 (female),
 * put its index to @param index
 * @return number of missing elements
 */
int findMissingSex(const std::vector<int>& sex, std::vector<int>* index);

/**
 * Remove i th element from @param val where i is stored in @param index
 * @return number of elements removed
 */
int removeByRowIndex(const std::vector<int>& index, Matrix* val);

/**
 * append a column @param val to the right of @param mat,
 * and set its label as @param label
 * @return 0 if success
 */
int appendToMatrix(const std::string& label, const std::vector<int> val,
                   Matrix* mat);

#endif /* _DATALOADER_H_ */
