#ifndef _LINEARREGRESSIONVT_H_
#define _LINEARREGRESSIONVT_H_

class Matrix;
class Vector;

class LinearRegressionVT {
 private:
  class LinearVTImpl;
  LinearVTImpl* impl;

 public:
  LinearRegressionVT();
  ~LinearRegressionVT();

  bool FitNullModel(const Matrix& Xnull, const Vector& y);

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
  LinearRegressionVT(const LinearRegressionVT& s);
  LinearRegressionVT& operator=(const LinearRegressionVT& s);
};

#endif /* _LINEARREGRESSIONVT_H_ */
