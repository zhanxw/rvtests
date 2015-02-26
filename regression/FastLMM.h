#ifndef _FASTLMM_H_
#define _FASTLMM_H_

#include <vector>

class EigenMatrix;
class Matrix;
class Vector;

class FastLMM{
 public:
  enum Test {
    SCORE = 0,
    LRT
  };
  enum Model {
    MLE = 0,
    REML
  };
 public: // Make this Impl public to make optimization function easy to write
  class Impl;
  Impl* impl;
 public:
  FastLMM(Test test, Model model);
  ~FastLMM();

  // @return 0 when success
  int FitNullModel(Matrix& Xnull, Matrix& y,
                   const EigenMatrix& kinshipU, const EigenMatrix& kinshipS);
  int TestCovariate(Matrix& Xnull, Matrix& y, Matrix& Xcol,
                    const EigenMatrix& kinshipU, const EigenMatrix& kinshipS);
  int CalculateUandV(Matrix& Xnull, Matrix& Y, Matrix& Xcol,
                     const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                     Matrix* uMat, Matrix* vMat);

  // NOTE: need to fit null model fit before calling this function
  double GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS);
  // NOTE: need to fit null model fit before calling this function
  double GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS, Matrix& Xcol);
  // NOTE: need to fit null model fit before calling this function
  double FastGetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS, Matrix& Xcol);
  double FastGetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS, Matrix& Xcol, int col);
  double GetPvalue();
  // for LRT Test
  double GetNullLogLikelihood() const;
  double GetAltLogLikelihood() const;
  // for Score Test
  double GetUStat() const;
  double GetVStat() const;
  int GetUMatrix(Matrix* u) const;
  int GetVMatrix(Matrix* v) const;
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
  // get covariance
  void GetCovZY(Matrix* zy);  // Cov(YZ) = Y' Simga^{-1} Z
  void GetCovZZ(Matrix* zz);  // Cov(ZZ) = Z' Simga^{-1} Z
  // transform genotype (e.g. in MetaCov)
  // x <- U' * ( x - center(x) )
  int TransformCentered(std::vector<double>* x,
                        const EigenMatrix& kinshipU,
                        const EigenMatrix& kinshipS);
  // @param out = sigma2_g * (lambda + delta) = sigma2_g * lambda + sigma2_e;
  int GetWeight(Vector* out) const;
  
};


#endif /* _FASTLMM_H_ */
