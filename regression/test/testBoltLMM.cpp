#include "base/IO.h"

#include "libVcf/PlinkOutputFile.h"
#include "libsrc/MathMatrix.h"
#include "libsrc/MathVector.h"
#include "regression/BoltLMM.h"
#include "regression/test/MatrixIO.h"

#include "base/SimpleMatrix.h"
#include "base/SimpleTimer.h"
#include "regression/EigenMatrix.h"
#include "regression/EigenMatrixInterface.h"

void normalize(Eigen::MatrixXf* mat);
void assignColumn(Matrix& G, int colIdx, Matrix* Xcol);
void writeToPlink(const Eigen::MatrixXf& y, const Eigen::MatrixXf& g,
                  const Eigen::MatrixXf& cov, const std::string& fn);

int main(int argc, char* argv[]) {
  fprintf(stderr, "Set omp threads = 12\n");
  omp_set_num_threads(12);
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

  LoadMatrix("input.bolt.g", G);
  LoadMatrix("input.bolt.y", Y);
  LoadMatrix("input.bolt.cov", Cov);

  BoltLMM bolt;

  EigenMatrix g;
  Eigen::MatrixXf y;
  Eigen::MatrixXf cov;

  G_to_Eigen(G, &g.mat);
  G_to_Eigen(Y, &y);
  G_to_Eigen(Cov, &cov);

  std::string fn = "input.bolt";
  writeToPlink(y, g.mat, cov, fn);

  // normalize(&g.mat);
  if (bolt.FitNullModel(fn, NULL) != 0) {
    fprintf(stderr, "Fit null model failed!");
    exit(1);
  }

  int M = g.mat.cols();
  Matrix Xcol;
  FILE* fp;
  if (argc == 1) {
    fp = fopen("output.bolt", "w");
  } else {
    fp = fopen(argv[1], "w");
  }
  fprintf(fp, "af\tu\tsqrt_v\teffect\tpval\n");
  for (int i = 0; i < M; ++i) {
    assignColumn(G, i, &Xcol);
    if (bolt.TestCovariate(Xcol)) {
      fprintf(stderr, "Test covariate failed!");
      exit(1);
    }
    fprintf(fp, "%g\t%g\t%g\t%g\t%g\n", bolt.GetAF(), bolt.GetU(),
            sqrt(bolt.GetV()), bolt.GetEffect(), bolt.GetPvalue());
  }
  fclose(fp);

#if 0      
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
    MultipleTraitLinearRegressionScoreTest mt(1024);

    bool ret = mt.FitNullModel(Cov, Y, tests);
    if (ret == false) {
      printf("Fit null model failed!\n");
      exit(1);
    }

    // for (int i = 0; i < N; ++i) {
    ret = mt.AddCovariate(G);
    if (ret == false) {
      printf("Add covariate failed!\n");
      exit(1);
    }
    ret = mt.TestCovariateBlock();
    if (ret == false) {
      printf("Test covariate block failed!\n");
      exit(1);
    }
    //}
    const Vector& pval = mt.GetPvalue(0);
    Print(pval);
    printf("\n");
  }

  // printf("T = %d\tN = %d\telapsed %.5f\n", T, N, t.stop());
#endif

  return 0;
}

void normalize(Eigen::MatrixXf* mat) {
  Eigen::MatrixXf& m = *mat;
  const int nrow = m.rows();
  const int ncol = m.cols();

  double sum = 0;
  int nonMiss = 0;
  double sum2 = 0;
  double avg, sdInv;
  for (int i = 0; i < ncol; i++) {
    // process column by column
    sum = 0;
    sum2 = 0;
    nonMiss = 0;
    for (int j = 0; j < nrow; ++j) {
      if (m(j, i) < 0) {
        continue;
      }
      sum += m(j, i);
      sum2 += m(j, i) * m(j, i);
      nonMiss++;
    }

    if (nonMiss == 0) {
      avg = 0.0;
      sdInv = 1.0;
    } else {
      avg = sum / nonMiss;
      sdInv = 1.0 / sqrt(sum2 / nonMiss - avg * avg);
    }

    for (int j = 0; j < nrow; ++j) {
      if (m(j, i) < 0) {
        m(j, i) = 0.;
      } else {
        m(j, i) = (m(j, i) - avg) * sdInv;
      }
    }
  }
}

void assignColumn(Matrix& G, int colIdx, Matrix* Xcol) {
  Matrix& x = *Xcol;
  x.Dimension(G.rows, 1);
  for (int i = 0; i < G.rows; ++i) {
    x[i][0] = G[i][colIdx];
  }
}

void writeToPlink(const Eigen::MatrixXf& y, const Eigen::MatrixXf& g,
                  const Eigen::MatrixXf& cov, const std::string& fn) {
  const int N = y.rows();
  const int M = g.cols();
  const int C = cov.cols();

  // say the fid/iid have the format: fid1/iid1, fid2/iid2, ...
  std::vector<std::string> fid;
  std::vector<std::string> iid;
  std::vector<double> pheno;
  char buf[128];
  for (int i = 0; i < N; ++i) {
    sprintf(buf, "fid%d", i);
    fid.push_back(buf);
    sprintf(buf, "iid%d", i);
    iid.push_back(buf);
    pheno.push_back(y(i, 0));
  }
  PlinkOutputFile pout(fn);
  // write .fam file
  pout.writeFAM(fid, iid, pheno);
  // write .bim file
  for (int i = 0; i < M; ++i) {
    sprintf(buf, "snp%d", i + 1);
    pout.writeBIM("1", buf, 0, i + 1, "A", "C");
  }
  // write .bed file
  SimpleMatrix m(M, N);  // marker by people
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      m[j][i] = g(i, j);
    }
  }
  pout.writeBED(&m, N, M);
  // finish
  pout.close();

  // write covar file
  std::string covFn = fn;
  covFn += ".covar";
  FILE* fp = fopen(covFn.c_str(), "wt");
  fprintf(fp, "FID\tIID");
  for (int i = 0; i < C; ++i) {
    sprintf(buf, "\tC%d", i + 1);
    fputs(buf, fp);
  }
  fputs("\n", fp);
  for (int i = 0; i < N; ++i) {
    fputs(fid[i].c_str(), fp);
    fputs("\t", fp);
    fputs(iid[i].c_str(), fp);
    for (int j = 0; j < C; ++j) {
      fputs("\t", fp);
      fprintf(fp, "%g", cov(i, j));
    }
    fputs("\n", fp);
  }
  fclose(fp);
}
