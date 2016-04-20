#include "MultipleTraitLinearRegressionScoreTest.h"

#include <map>
#include <string>
#include <vector>

#include "gsl/gsl_cdf.h"  // use gsl_cdf_chisq_Q

#include "Eigen/Dense"
#include "Eigen/Cholesky"  // ldlt

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
  EMatVec Uyz;
  EMatVec Ugz;
  EMatVec Uyg;
  std::vector<bool> hasCovariate;
  std::vector<std::vector<bool> >
      missingIndex;  // store whether the elements of Y[i] or Z[i, .] is missing
  std::vector<double> sigma2;
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
  (*m).colwise().normalize();
  // rowwise() /= (*m).colwise().norm().array();
}

MultipleTraitLinearRegressionScoreTest::
    MultipleTraitLinearRegressionScoreTest() {
  this->work = new MultipleTraitLinearRegressionScoreTestInternal;
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
  w.hasCovariate.resize(tests.size());
  w.missingIndex.resize(tests.size());
  w.Uyz.resize(tests.size());
  w.Ugz.resize(tests.size());
  w.Uyg.resize(tests.size());
  w.sigma2.resize(tests.size());
  w.nTest = tests.size();
  ustat.Dimension(tests.size());
  vstat.Dimension(tests.size());
  pvalue.Dimension(tests.size());

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

  // NOTE: need to handle (1) no covariates (2) after removing missing values,
  // there is no value to test
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
      w.ZZinv[i].noalias() = (w.Z[i].transpose() * w.Z[i]).ldlt().solve(
          EMat::Identity(w.Z[i].cols(), w.Z[i].cols()));
      w.Uyz[i].noalias() = w.Z[i].transpose() * w.Y[i];
      w.sigma2[i] = (w.Y[i].transpose() * w.Y[i] -
                     w.Uyz[i].transpose() * w.ZZinv[i] * w.Uyz[i])(0, 0) /
                    w.Y[i].rows();
    } else {
      w.sigma2[i] = w.Y[i].col(0).squaredNorm() / w.Y[i].rows();
    }
  }

  return true;
}
bool MultipleTraitLinearRegressionScoreTest::TestCovariate(Matrix& g) {
  MultipleTraitLinearRegressionScoreTestInternal& w = *this->work;

  for (int i = 0; i < w.nTest; ++i) {
    // Convert g to suitable g matrix
    const int nSample = w.Y[i].rows();
    w.G[i].resize(nSample, 1);
    int idx = 0;
    for (int j = 0; j < (int)w.missingIndex[i].size(); ++j) {
      if (!w.missingIndex[i][j]) {
        w.G[i](idx, 0) = g[j][0];
        ++idx;
      }
    }

    // center and scale g
    scale(&w.G[i]);

    // calculate Ugz, Uyg
    if (w.hasCovariate[i]) {
      w.Ugz[i].noalias() = w.Z[i].transpose() * w.G[i];  // C by 1
    }
    w.Uyg[i].noalias() = w.G[i].transpose() * w.Y[i];

    // calculate Ustat, Vstat
    if (w.hasCovariate[i]) {
      ustat[i] =
          (w.Uyg[i] - w.Ugz[i].transpose() * w.ZZinv[i] * w.Uyz[i])(0, 0);
      vstat[i] = (w.G[i].transpose() * w.G[i] -
                  w.Ugz[i].transpose() * w.ZZinv[i] * w.Ugz[i])(0, 0);
    } else {  // no covariate
      ustat[i] = w.Uyg[i](0, 0);
      vstat[i] =
          w.G[i].col(0).squaredNorm();  // (w.G[i].transpose() * w.G[i])(0, 0);
    }
    vstat[i] *= w.sigma2[i];

    // calculat p-value
    if (vstat[i] == 0.) {
      pvalue[i] = NAN;
    } else {
      double stat = ustat[i] * ustat[i] / vstat[i];
      pvalue[i] = gsl_cdf_chisq_Q(stat, 1.0);
    }
  }
  return true;
}
