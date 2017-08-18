#ifndef _BOLTLMM_H_
#define _BOLTLMM_H_

#include <string>
#include <vector>

class Matrix;
class FloatMatrixRef;

class BoltLMM {
  class BoltLMMImpl;

 public:
  BoltLMM();
  ~BoltLMM();
  /**
   * Null model will estimate the heritability
   * @return 0 when success
   * @param prefix   PLINK file prefix, also, "${prefix}.covar" will store
   * @param pheno    provide phenotypes; NULL => read phentoypes from PLINK file
   * covaraites and covariate file usually does not have intercept
   */
  int FitNullModel(const std::string& prefix, Matrix* pheno);
  // test @param Xcol, a column matrix
  int TestCovariate(Matrix& Xcol);

  // get results
  double GetAF();
  double GetU();
  double GetV();
  double GetEffect();
  double GetPvalue();

  // calculate covariances
  void GetCovXX(const std::vector<double>& g1, const std::vector<double>& g2,
                double* out);
  void GetCovXX(const FloatMatrixRef& g1, const FloatMatrixRef& g2, float* out);

 private:
  // don't copy
  BoltLMM(const BoltLMM&);
  BoltLMM& operator=(const BoltLMM&);

 private:
  BoltLMMImpl* impl_;
};

#endif /* _BOLTLMM_H_ */
