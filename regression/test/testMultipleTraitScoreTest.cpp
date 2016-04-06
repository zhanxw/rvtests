#include "IO.h"

#include "MathMatrix.h"
#include "MathVector.h"
#include "MultipleTraitLinearRegressionScoreTest.h"
#include "MatrixIO.h"

#include "Formula.h"
#include "SimpleTimer.h"

int main(int argc, char* argv[]) {
  // int T = 10;
  // int N = 10000;
  // if (argc >= 2) {
  //   T = atoi(argv[1]);
  // }
  // if (argc >= 3) {
  //   N = atoi(argv[2]);
  // }

  Matrix G;
  Matrix Y;
  Matrix Cov;

  LoadMatrix("input.mt.g", G);
  LoadMatrix("input.mt.y", Y);
  LoadMatrix("input.mt.cov", Cov);
  Cov.SetColumnLabel(0, "c1");
  Cov.SetColumnLabel(1, "c2");
  Y.SetColumnLabel(0, "y1");
  Y.SetColumnLabel(1, "y2");
  Y.SetColumnLabel(2, "y3");

  FormulaVector tests;
  // std::vector<std::vector<std::string> > v(2);
  {
    const char* tp1[] = {"y1"};
    const char* tc1[] = {"c1"};
    std::vector<std::string> p1(tp1, tp1 + 1);
    std::vector<std::string> c1(tc1, tc1 + 1);
    // v[0] = p1;
    // v[1] = c1;
    // tests.push_back(v);
    tests.add(p1, c1);
  }

  {
    const char* tp1[] = {"y2"};
    const char* tc1[] = {"c2"};
    std::vector<std::string> p1(tp1, tp1 + 1);
    std::vector<std::string> c1(tc1, tc1 + 1);
    // v[0] = p1;
    // v[1] = c1;
    // tests.push_back(v);
    tests.add(p1, c1);
  }

  // for (int i = 0; i < T - 2; ++i) {
  {
    const char* tp1[] = {"y2"};
    const char* tc1[] = {"c1", "c2"};
    std::vector<std::string> p1(tp1, tp1 + 1);
    std::vector<std::string> c1(tc1, tc1 + 2);
    // v[0] = p1;
    // V[1] = c1;
    // tests.push_back(v);
    tests.add(p1, c1);
  }

  AccurateTimer t;
  {
    MultipleTraitLinearRegressionScoreTest mt;

    bool ret = mt.FitNullModel(Cov, Y, tests);
    if (ret == false) {
      printf("Fit null model failed!\n");
      exit(1);
    }

    // for (int i = 0; i < N; ++i) {
    ret = mt.TestCovariate(G);
    if (ret == false) {
      printf("Test covariate failed!\n");
      exit(1);
    }
    //}
    Vector& pval = mt.GetPvalue();
    Print(pval);
    printf("\n");
  }

  // printf("T = %d\tN = %d\telapsed %.5f\n", T, N, t.stop());

  return 0;
}
