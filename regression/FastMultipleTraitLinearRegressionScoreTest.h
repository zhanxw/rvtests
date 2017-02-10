#ifndef FASTMULTIPLETRAITSCORETEST_H
#define FASTMULTIPLETRAITSCORETEST_H

#include "libsrc/MathMatrix.h"

class FastMultipleTraitLinearRegressionScoreTestInternal;
class FormulaVector;

class FastMultipleTraitLinearRegressionScoreTest {
 public:
  FastMultipleTraitLinearRegressionScoreTest(int blockSize);
  virtual ~FastMultipleTraitLinearRegressionScoreTest();
  bool FitNullModel(Matrix& cov, Matrix& y, const FormulaVector& tests);
  bool AddGenotype(const Matrix& g);
  bool TestCovariateBlock();
  const Vector& GetPvalue(int i) const { return this->pvalue[i]; };
  const Vector& GetU(int i) const { return this->ustat[i]; };
  const Vector& GetV(int i) const { return this->vstat[i]; };
  void flush() { resultLength = 0; };

 private:
  FastMultipleTraitLinearRegressionScoreTest(
      const FastMultipleTraitLinearRegressionScoreTest&);
  FastMultipleTraitLinearRegressionScoreTest& operator=(
      const FastMultipleTraitLinearRegressionScoreTest&);

 private:
  Matrix freq;    // [ resultLen x nTest ]
  Matrix ustat;   // [ resultLen x nTest ]
  Matrix vstat;   // [ resultLen x nTest ]
  Matrix pvalue;  // [ resultLen x nTest ]
  FastMultipleTraitLinearRegressionScoreTestInternal*
      work;          // store working data
  int blockSize;     // unit of grouped computational units
  int resultLength;  // how many results are available (always < blockSize)
};

#endif /* FASTMULTIPLETRAITSCORETEST_H */
