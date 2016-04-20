#ifndef MULTIPLETRAITSCORETEST_H
#define MULTIPLETRAITSCORETEST_H

#include <string>
#include <vector>
#include "MathMatrix.h"

class MultipleTraitLinearRegressionScoreTestInternal;
class FormulaVector;

class MultipleTraitLinearRegressionScoreTest {
 public:
  MultipleTraitLinearRegressionScoreTest();
  virtual ~MultipleTraitLinearRegressionScoreTest();
  bool FitNullModel(Matrix& cov, Matrix& y,
                    const FormulaVector & tests);
  bool TestCovariate(Matrix& g);
  Vector& GetPvalue() { return this->pvalue; };

 private:
  Vector ustat;
  Vector vstat;  
  Vector pvalue;
  MultipleTraitLinearRegressionScoreTestInternal* work;  // store working data
};

#endif /* MULTIPLETRAITSCORETEST_H */
