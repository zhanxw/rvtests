#ifndef _DATALOADER_H_
#define _DATALOADER_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#include "base/MathMatrix.h"
#include "base/MathVector.h"
#include "base/SimpleMatrix.h"
#include "regression/Formula.h"

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
  int addFittedParameter(const std::string& name, double beta, double seBeta);
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
  const SimpleMatrix& getEstimation() const {
    return this->fittedResidualModel;
  };

 private:
  // for multiple traits
  std::vector<SimpleMatrix> phenotypes;  // sample by traits
  std::vector<SimpleMatrix> covariates;  // sample by covariates
  FormulaVector formula;

  SimpleMatrix& phenotype;  // sample by traits
  SimpleMatrix& covariate;  // sample by covariates
  bool binaryPhenotype;
  std::vector<int> sex;  // plink coded genders

  SimpleMatrix fittedResidualModel;  // store estimates (beta, se(beta)) for
                                     // model y ~ cov,

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
};  // end DataLoader

#endif /* _DATALOADER_H_ */
