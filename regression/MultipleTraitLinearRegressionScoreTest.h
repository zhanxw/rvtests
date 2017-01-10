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
  bool AddGenotype(const Matrix& g);
  bool TestCovariateBlock();
  const Vector& GetPvalue(int i) const { return this->pvalue[i]; };
  const Vector& GetU(int i) const { return this->ustat[i]; };
  const Vector& GetV(int i) const { return this->vstat[i]; };
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
