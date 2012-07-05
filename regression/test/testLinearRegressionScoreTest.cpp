#include "IO.h"

#include "MathMatrix.h"
#include "MathVector.h"
#include "LinearRegressionScoreTest.h"
#include "MatrixIO.h"

int main(int argc, char *argv[])
{
  Vector Y;
  Matrix X;

  LoadVector("input.y", Y);
  LoadMatrix("input.x", X);

  {
    Matrix x;
    Vector y;
    x = X;
    y = Y;

    LinearRegressionScoreTest lrst;
    if ( lrst.FitLinearModel(x, y, 1) == false) {
      fprintf(stderr, "Fitting failed!\n");
      return -1;
    }

    fprintf(stdout, "score_p\t");
    double score_p =lrst.GetPvalue();
    Print(score_p);
    fputc('\n', stdout);
  }

  {
    Matrix intercept;
    Matrix x;
    extractColumn(X, 0, &intercept);
    extractColumn(X, 1, &x);
    Vector y;
    y = Y;
    
    LinearRegressionScoreTest lrst;
    bool ret = lrst.FitNullModel(intercept, y);
    assert(ret);
    ret = lrst.TestCovariate(intercept, y, x);
    assert(ret);

    fprintf(stdout, "score_p\t");
    double score_p =lrst.GetPvalue();
    Print(score_p);
    fputc('\n', stdout);
  }

  {
    Matrix intercept;
    Vector x;
    extractColumn(X, 0, &intercept);
    extractColumn(X, 1, &x);
    Vector y;
    y = Y;
    
    LinearRegressionScoreTest lrst;
    bool ret = lrst.FitNullModel(intercept, y);
    assert(ret);
    ret = lrst.TestCovariate(intercept, y, x);
    assert(ret);

    fprintf(stdout, "score_p\t");
    double score_p =lrst.GetPvalue();
    Print(score_p);
    fputc('\n', stdout);
  }

  {
    Vector x;
    extractColumn(X, 1, &x);
    Vector y;
    y = Y;
    
    LinearRegressionScoreTest lrst;
    bool ret = lrst.TestCovariate(x, y);
    assert(ret);

    fprintf(stdout, "score_p\t");
    double score_p =lrst.GetPvalue();
    Print(score_p);
    fputc('\n', stdout);
  }

  {
    Matrix x;
    extractColumn(X, 1, &x);
    Vector y;
    y = Y;
    
    LinearRegressionScoreTest lrst;
    bool ret = lrst.TestCovariate(x, y);
    assert(ret);

    fprintf(stdout, "score_p\t");
    double score_p =lrst.GetPvalue();
    Print(score_p);
    fputc('\n', stdout);
  }
  return 0;
};

