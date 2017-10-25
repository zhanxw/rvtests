#ifndef _METACOV_H_
#define _METACOV_H_

#include <vector>

class EigenMatrix;
class Matrix;
class Vector;

/**
 * This class is to help calculate covariance in family strucutre
 */
class MetaCov {
 public:  // Make this Impl public to make optimization function easy to write
  class Impl;
  Impl* impl;

 public:
  MetaCov();
  ~MetaCov();

  // @return 0 when success
  int FitNullModel(Matrix& Xnull, Matrix& y, const EigenMatrix& kinshipU,
                   const EigenMatrix& kinshipS);
  // U' * ( x - center(x) )
  int TransformCentered(std::vector<double>* x, const EigenMatrix& kinshipU,
                        const EigenMatrix& kinshipS);
  // U' * x
  int Transform(std::vector<double>* x, const EigenMatrix& kinshipU,
                const EigenMatrix& kinshipS);

  // @param out = sigma2_g * (lambda + delta) = sigma2_g * lambda + sigma2_e;
  int GetWeight(Vector* out);

  // @parm out = G' (\Sigma^{-1} - ...) * G
  //                     ^ scaled sigma
  // NOTE: here assume @param g1 and @param g2 is transformed. e.g. U' * g
  void GetCovXX(const std::vector<double>& g1, const std::vector<double>& g2,
                const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                double* out);

  // @parm out = G' \Sigma^{-1} Z, @param is a row vector
  // NOTE: here assume @param g is transformed. e.g. U' * g
  void GetCovXZ(const std::vector<double>& g, const EigenMatrix& kinshipU,
                const EigenMatrix& kinshipS, std::vector<double>* out);

  // @parm out = Z' \Sigma^{-1} Z
  void GetCovZZ(Matrix* out);
};

#endif /* _METACOV_H_ */
