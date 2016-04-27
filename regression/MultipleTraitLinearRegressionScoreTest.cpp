#include "MultipleTraitLinearRegressionScoreTest.h"

#include <map>
#include <string>
#include <vector>

#include "gsl/gsl_cdf.h"  // use gsl_cdf_chisq_Q

#include "Eigen/Cholesky"  // ldlt
#include "Eigen/Dense"

#include "Formula.h"

typedef Eigen::MatrixXf EMat;
typedef std::vector<EMat> EMatVec;

class MultipleTraitLinearRegressionScoreTestInternal {
 public:
  // these are counts before removing missing items
  int N;  // sample
  int T;  // trait
  int M;  // marker
  int C;  // covariates
  int nTest;

  // following items are all of length = nTests
  EMatVec Y;
  EMatVec Z;
  EMatVec G;
  EMatVec ZZinv;
  EMatVec L;
  EMatVec Uyz;
  EMatVec Ugz;
  EMatVec Uyg;
  std::vector<bool> hasCovariate;
  std::vector<std::vector<bool> >
      missingIndex;  // store whether the elements of Y[i] or Z[i, .] is missing
  std::vector<double> sigma2;
  EMat ustat;
  EMat vstat;
};

void makeColNameToDict(Matrix& m, std::map<std::string, int>* dict) {
  std::map<std::string, int>& d = *dict;
  for (int i = 0; i < m.cols; ++i) {
    d[m.GetColumnLabel(i)] = i;
  }
}

void makeMatrix(Matrix& m, const std::vector<int>& index, EMat* out) {
  (*out).resize(m.rows, index.size());
  for (int i = 0; i < m.rows; ++i) {
    for (size_t j = 0; j < index.size(); ++j) {
      const int idx = index[j];
      (*out)(i, j) = m[i][idx];
    }
  }
}

bool hasMissingInRow(const EMat& m, int r) {
  const int nc = m.cols();
  for (int i = 0; i < nc; ++i) {
    if (isnan(m(r, i))) {
      return true;
    }
  }
  return false;
}

/**
 * Remove i th row in @param m if i th element in @param missingIndicator is
 * true
 */
void removeRow(const std::vector<bool>& missingIndicator, EMat* m) {
  const int nr = (*m).rows();
  int idx = 0;
  for (int i = 0; i < nr; ++i) {
    if (missingIndicator[i]) continue;
    if (idx != i) {
      (*m).row(idx) = (*m).row(i);
    }
    ++idx;
  }
  (*m).conservativeResize(idx, (*m).cols());
}

void scale(EMat* m) {
  (*m).rowwise() -= (*m).colwise().sum() / (*m).rows();
  // (*m).colwise().normalize();
}

MultipleTraitLinearRegressionScoreTest::MultipleTraitLinearRegressionScoreTest(
    int blockSize) {
  this->work = new MultipleTraitLinearRegressionScoreTestInternal;
  this->blockSize = blockSize;
  this->resultLength = 0;
}

MultipleTraitLinearRegressionScoreTest::
    ~MultipleTraitLinearRegressionScoreTest() {
  if (this->work) {
    delete this->work;
    this->work = NULL;
  }
}
bool MultipleTraitLinearRegressionScoreTest::FitNullModel(
    Matrix& cov, Matrix& pheno, const FormulaVector& tests) {
  MultipleTraitLinearRegressionScoreTestInternal& w = *this->work;
  // set some values
  w.N = cov.rows;
  w.T = pheno.cols;
  w.C = cov.cols;
  w.M = -1;

  w.Y.resize(tests.size());
  w.Z.resize(tests.size());
  w.G.resize(tests.size());
  w.ZZinv.resize(tests.size());
  w.L.resize(tests.size());
  w.hasCovariate.resize(tests.size());
  w.missingIndex.resize(tests.size());
  w.Uyz.resize(tests.size());
  w.Ugz.resize(tests.size());
  w.Uyg.resize(tests.size());
  w.sigma2.resize(tests.size());
  w.nTest = tests.size();
  w.ustat.resize(blockSize, tests.size());
  w.vstat.resize(blockSize, tests.size());
  ustat.Dimension(blockSize, tests.size());
  vstat.Dimension(blockSize, tests.size());
  pvalue.Dimension(blockSize, tests.size());

  // create dict (key: phenotype/cov name, val: index)
  std::map<std::string, int> phenoDict;
  std::map<std::string, int> covDict;
  makeColNameToDict(pheno, &phenoDict);
  makeColNameToDict(cov, &covDict);

  // create Y, Z
  std::vector<std::string> phenoName;
  std::vector<std::string> covName;
  std::vector<int> phenoCol;
  std::vector<int> covCol;

  // NOTE: need to handle two corner cases:
  // (1) no covariates
  // (2) after removing missing values, there is no value to test
  for (int i = 0; i < w.nTest; ++i) {
    phenoName = tests.getPhenotype(i);
    phenoCol.clear();
    phenoCol.push_back(phenoDict[phenoName[0]]);
    covName = tests.getCovariate(i);
    covCol.clear();
    for (size_t j = 0; j != covName.size(); ++j) {
      if (covName[j] == "1") {
        continue;
      }
      assert(covDict.count(covName[j]));
      covCol.push_back(covDict[covName[j]]);
    }

    w.hasCovariate[i] = covCol.size() > 0;
    makeMatrix(pheno, phenoCol, &w.Y[i]);
    if (w.hasCovariate[i]) {
      makeMatrix(cov, covCol, &w.Z[i]);
    }

    // create index to indicate missingness
    w.missingIndex[i].resize(w.N);
    for (int j = 0; j < w.N; ++j) {
      if (hasMissingInRow(w.Y[i], j)) {
        w.missingIndex[i][j] = true;
        continue;
      } else {
        if (w.hasCovariate[i] && hasMissingInRow(w.Z[i], j)) {
          w.missingIndex[i][j] = true;
          continue;
        }
      }
      w.missingIndex[i][j] = false;
    }
    removeRow(w.missingIndex[i], &w.Y[i]);
    removeRow(w.missingIndex[i], &w.Z[i]);

    // center and scale Y, Z
    scale(&w.Y[i]);
    scale(&w.Z[i]);

    // calcualte Uzy, inv(Z'Z)
    if (w.hasCovariate[i]) {
      w.ZZinv[i].noalias() =
          (w.Z[i].transpose() * w.Z[i])
              .ldlt()
              .solve(EMat::Identity(w.Z[i].cols(), w.Z[i].cols()));
      Eigen::LLT<Eigen::MatrixXf> lltOfA(w.ZZinv[i]);
      // L * L' = A
      w.L[i] = lltOfA.matrixL();
      w.Uyz[i].noalias() = w.Z[i].transpose() * w.Y[i];
      w.sigma2[i] = (w.Y[i].transpose() * w.Y[i] -
                     w.Uyz[i].transpose() * w.ZZinv[i] * w.Uyz[i])(0, 0) /
                    w.Y[i].rows();
    } else {
      w.sigma2[i] = w.Y[i].col(0).squaredNorm() / w.Y[i].rows();
    }

    // initialize G
    const int nSample = w.Y[i].rows();
    if (w.G[i].cols() != blockSize) {
      w.G[i].resize(nSample, blockSize);
    }
  }

  return true;
}

bool MultipleTraitLinearRegressionScoreTest::AddCovariate(const Matrix& g) {
  MultipleTraitLinearRegressionScoreTestInternal& w = *this->work;
  assert(resultLength < blockSize);
  for (int i = 0; i < w.nTest; ++i) {
    // Convert g to suitable g matrix
    int idx = 0;
    const std::vector<bool>& missingIndex = w.missingIndex[i];
    const int n = missingIndex.size();
    EMat& G = w.G[i];

    for (int j = 0; j < n; ++j) {
      if (!missingIndex[j]) {
        G(idx, resultLength) = g[j][0];
        ++idx;
      }
    }
  }
  resultLength++;
  return true;
}

bool MultipleTraitLinearRegressionScoreTest::TestCovariateBlock() {
  MultipleTraitLinearRegressionScoreTestInternal& w = *this->work;
  for (int i = 0; i < w.nTest; ++i) {
    // delcare const variables
    EMat& G = w.G[i];
    const EMat& Z = w.Z[i];
    const EMat& Y = w.Y[i];
    EMat& Ugz = w.Ugz[i];
    EMat& Uyg = w.Uyg[i];
    const EMat& Uyz = w.Uyz[i];
    const bool& hasCovariate = w.hasCovariate[i];

    // center and scale g
    scale(&G);

    // calculate Ugz, Uyg
    if (hasCovariate) {
      Ugz.noalias() = Z.transpose() * G;  // C by blockSize
    }
    Uyg.noalias() = G.transpose() * Y;  // blockSize by T=1

    // calculate Ustat, Vstat
    if (hasCovariate) {
      w.ustat.col(i).noalias() =
          (Uyg - Ugz.transpose() * w.ZZinv[i] * Uyz);  // blockSize by T=1
      w.vstat.col(i).noalias() =
          (G.array().square() - (w.L[i].transpose() * Ugz).array().square())
              .matrix()
              .colwise()
              .sum()
              .transpose();  // blockSize by 1
    } else {                 // no covariate
      w.ustat.col(i).noalias() = Uyg;
      w.vstat.col(i).noalias() =
          G.array().square().matrix().colwise().sum().transpose();  // blockSize
                                                                    // by 1
    }
    w.vstat.col(i) *= w.sigma2[i];
  }
  // assign and calculat p-value
  for (int j = 0; j < blockSize; ++j) {
    for (int i = 0; i < w.nTest; ++i) {
      this->ustat[j][i] = w.ustat(j, i);
      this->vstat[j][i] = w.vstat(j, i);

      if (w.vstat(j, i) == 0.) {
        pvalue[j][i] = NAN;
      } else {
        double stat = ustat[j][i] * ustat[j][i] / vstat[j][i];
        pvalue[j][i] = gsl_cdf_chisq_Q(stat, 1.0);
      }
    }
  }  // end for i
  return true;
}
