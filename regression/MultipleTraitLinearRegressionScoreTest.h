#ifndef MULTIPLETRAITSCORETEST_H
#define MULTIPLETRAITSCORETEST_H

#include <string>
#include <vector>
#include "MathMatrix.h"

class MultipleTraitLinearRegressionScoreTestInternal;
class FormulaVector;

class MultipleTraitLinearRegressionScoreTest {
 public:
  MultipleTraitLinearRegressionScoreTest(int blockSize);
  virtual ~MultipleTraitLinearRegressionScoreTest();
  bool FitNullModel(Matrix& cov, Matrix& y, const FormulaVector& tests);
  bool AddCovariate(Matrix& g);
  bool TestCovariateBlock();
  Vector& GetPvalue(int i) { return this->pvalue[i]; };
  Vector& GetU(int i) { return this->ustat[i]; };
  Vector& GetV(int i) { return this->vstat[i]; };
  void flush() { resultLength = 0; };

 private:
  MultipleTraitLinearRegressionScoreTest(
      const MultipleTraitLinearRegressionScoreTest&);
  MultipleTraitLinearRegressionScoreTest& operator=(
      const MultipleTraitLinearRegressionScoreTest&);

 private:
  Matrix ustat;
  Matrix vstat;
  Matrix pvalue;
  MultipleTraitLinearRegressionScoreTestInternal* work;  // store working data
  int blockSize;     // unit of grouped computational units
  int resultLength;  // how many results are available
};

#endif /* MULTIPLETRAITSCORETEST_H */
