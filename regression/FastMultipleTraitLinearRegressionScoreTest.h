#ifndef FASTMULTIPLETRAITSCORETEST_H
#define FASTMULTIPLETRAITSCORETEST_H

#include "base/MathMatrix.h"

class FastMultipleTraitLinearRegressionScoreTestInternal;
class FormulaVector;

class FastMultipleTraitLinearRegressionScoreTest {
 public:
  FastMultipleTraitLinearRegressionScoreTest(int blockSize);
  virtual ~FastMultipleTraitLinearRegressionScoreTest();
  bool FitNullModel(const Matrix& cov, const Matrix& y,
                    const FormulaVector& tests);
  bool AddGenotype(const Matrix& g);
  bool TestCovariateBlock();
  const Matrix& GetPvalue() const { return this->pvalue; };
  const Matrix& GetU() const { return this->ustat; };
  const Matrix& GetV() const { return this->vstat; };
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
