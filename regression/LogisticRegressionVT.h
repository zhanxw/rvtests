#ifndef _LOGISTICREGRESSIONVT_H_
#define _LOGISTICREGRESSIONVT_H_

// #include "MathMatrix.h"

class Matrix;
class Vector;

class LogisticRegressionVT {
 private:
  class LogisticVTImpl;
  LogisticVTImpl* impl;

 public:
  LogisticRegressionVT();
  ~LogisticRegressionVT();

  bool FitNullModel(const Matrix& Xnull, const Vector& y, int nRound = 100);

  /**
   * Test H0: \beta = 0  (\beta is multiple dimension).
   * y ~ \beta * Xcol + \gamma * Xnull
   * @return false if not working
   */
  bool TestCovariate(const Matrix& Xnull, const Vector& y, const Matrix& Xcol);

  int GetIndexMax();  // return index to the maximum t
  Matrix& GetU();
  Matrix& GetV();
  Matrix& GetT();
  Matrix& GetCov();  // return cov(U)
  double GetPvalue() const;
  double GetEffect(int index) const;

 private:
  // don't copy
  LogisticRegressionVT(const LogisticRegressionVT& s);
  LogisticRegressionVT& operator=(const LogisticRegressionVT& s);
};

#endif /* _LOGISTICREGRESSIONVT_H_ */
