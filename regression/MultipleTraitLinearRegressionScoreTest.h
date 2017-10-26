#ifndef MULTIPLETRAITSCORETEST_H
#define MULTIPLETRAITSCORETEST_H

#include <string>
#include <vector>
#include "base/MathMatrix.h"

class MultipleTraitLinearRegressionScoreTestInternal;
class FormulaVector;

class MultipleTraitLinearRegressionScoreTest {
 public:
  MultipleTraitLinearRegressionScoreTest(int blockSize);
  virtual ~MultipleTraitLinearRegressionScoreTest();
  bool FitNullModel(Matrix& cov, Matrix& y, const FormulaVector& tests);
  bool AddGenotype(const Matrix& g);
  bool TestCovariateBlock();
  const Matrix& GetPvalue() const { return this->pvalue; };
  const Matrix& GetU() const { return this->ustat; };
  const Matrix& GetV() const { return this->vstat; };
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
  // store grouping results
  std::vector<std::vector<int> > group;  // record grouping results of tests
  int groupSize;
};

#endif /* MULTIPLETRAITSCORETEST_H */
