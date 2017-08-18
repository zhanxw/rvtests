#ifndef _FASTLMM_H_
#define _FASTLMM_H_

#include <vector>
#include "regression/MatrixRef.h"

class EigenMatrix;
class Matrix;
class Vector;

class FastLMM {
 public:
  enum Test { SCORE = 0, LRT };
  enum Model { MLE = 0, REML };

 public:  // Make this Impl public to make optimization function easy to write
  class Impl;
  Impl* impl;

 public:
  FastLMM(Test test, Model model);
  ~FastLMM();

  // @return 0 when success
  int FitNullModel(Matrix& Xnull, Matrix& y, const EigenMatrix& kinshipU,
                   const EigenMatrix& kinshipS);
  // test @param Xcol
  int TestCovariate(Matrix& Xnull, Matrix& y, Matrix& Xcol,
                    const EigenMatrix& kinshipU, const EigenMatrix& kinshipS);
  // U = (X_centered)' * Inverse(Sigma) * Y_res
  // V = (X_centered)' * Inverse(Sigma_ X|Y) * X_centered
  int CalculateUandV(Matrix& Xnull, Matrix& Y, Matrix& Xcol,
                     const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                     Matrix* uMat, Matrix* vMat);

  // NOTE: need to fit null model fit before calling this function
  double GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS);
  // NOTE: need to fit null model fit before calling this function
  double GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
               Matrix& Xcol);
  // NOTE: need to fit null model fit before calling this function
  double FastGetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                   Matrix& Xcol);
  double FastGetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                   Matrix& Xcol, int col);
  double GetPvalue();
  // for LRT Test
  double GetNullLogLikelihood() const;
  double GetAltLogLikelihood() const;
  // for Score Test
  double GetUStat() const;    // X' \Simga^{-1} Y
  double GetVStat() const;    // X' \Simga^{-1} X
  double GetEffect() const;   // U/V
  double GetSE() const;       // 1/sqrt(V)
  double GetSigmaG2() const;  // sigma_g^2
  double GetSigmaE2() const;  // sigma_e^2
  double GetDelta() const;    // delta = sigma2_e / sigma2_g
  void GetBeta(EigenMatrix* beta) const;
  // get estimates from null
  void GetNullCovEst(Vector* beta);
  void GetNullCovB(Matrix* betaCov);
  // get scaled factors
  double GetSigmaK();
  double GetSigma1();
  // get covariance. Z: covariates, Y: response, X: genotype
  void GetCovZY(Matrix* zy);  // Cov(YZ) = Y' Simga^{-1} Z
  void GetCovZZ(Matrix* zz);  // Cov(ZZ) = Z' Simga^{-1} Z
  // NOTE: here assume @param g1 and @param g2 is transformed. e.g. U' * g
  void GetCovXX(const std::vector<double>& g1, const std::vector<double>& g2,
                const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                double* out);
  void GetCovXX(FloatMatrixRef& g1, FloatMatrixRef& g2,
                const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                float* out);
  // NOTE: here assume @param g is transformed. e.g. U' * g
  void GetCovXZ(const std::vector<double>& g, const EigenMatrix& kinshipU,
                const EigenMatrix& kinshipS, std::vector<double>* out);
  void GetCovXZ(FloatMatrixRef& g, const EigenMatrix& kinshipU,
                const EigenMatrix& kinshipS, FloatMatrixRef& out);
  // transform genotype (e.g. in MetaCov)
  // x <- U' * ( x - center(x) )
  int TransformCentered(std::vector<double>* x, const EigenMatrix& kinshipU,
                        const EigenMatrix& kinshipS);
  int TransformCentered(FloatMatrixRef& x, const EigenMatrix& kinshipU,
                        const EigenMatrix& kinshipS);
  // transform genotype (e.g. in MetaCov)
  // x <- U' * x
  int Transform(std::vector<double>* x, const EigenMatrix& kinshipU,
                const EigenMatrix& kinshipS);
  int Transform(FloatMatrixRef& geno, const EigenMatrix& kinshipU,
                const EigenMatrix& kinshipS);

  // @param out = sigma2_g * (lambda + delta) = sigma2_g * lambda + sigma2_e;
  int GetWeight(Vector* out) const;

  void disableCenterGenotype();
};

#endif /* _FASTLMM_H_ */
