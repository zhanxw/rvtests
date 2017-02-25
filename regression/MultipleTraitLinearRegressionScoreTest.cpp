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
  EMatVec Uyz;
  EMatVec Ugz;
  EMatVec Uyg;
  std::vector<bool> hasCovariate;
  std::vector<std::vector<bool> >
      missingIndex;  // store whether the elements of Y[i] or Z[i, .] is missing
  std::vector<double> sigma2;
  EMatVec ustat;
  EMatVec vstat;
  // grouped values
  EMatVec groupedY;
  EMatVec groupedZ;
  std::vector<bool> groupedHasCovariate;
  EMatVec groupedUyz;
  EMatVec groupedZZinv;
  EMatVec groupedL;
};

/// Column names of @param m are stored in @param dict
void makeColNameToDict(Matrix& m, std::map<std::string, int>* dict) {
  std::map<std::string, int>& d = *dict;
  for (int i = 0; i < m.cols; ++i) {
    d[m.GetColumnLabel(i)] = i;
  }
}

/// Extract columns, specified by @param index, from @param m, to @param out
void makeMatrix(Matrix& m, const std::vector<int>& index, EMat* out) {
  (*out).resize(m.rows, index.size());
  for (int i = 0; i < m.rows; ++i) {
    for (size_t j = 0; j < index.size(); ++j) {
      const int idx = index[j];
      (*out)(i, j) = m[i][idx];
    }
  }
}

/// Check whether row @param r has missing elements (denoted by NAN).
bool hasMissingInRow(const EMat& m, int r) {
  const int nc = m.cols();
  for (int i = 0; i < nc; ++i) {
    if (std::isnan(m(r, i))) {
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

/**
 * Scale matrix @param m to have zero mean, and unit variance.
 * NOTE: assume @param m does not have missing values
 */
void scale(EMat* m) {
  (*m).rowwise() -= (*m).colwise().sum() / (*m).rows();
  // (*m).colwise().normalize();
}

/**
 * Get a string of 0 and 1, it represents @param v
 */
std::string toString(const std::vector<bool> v) {
  std::string s;
  for (size_t i = 0; i != v.size(); ++i) {
    if (v[i]) {
      s.push_back('1');
    } else {
      s.push_back('0');
    }
  }
  return s;
}

MultipleTraitLinearRegressionScoreTest::MultipleTraitLinearRegressionScoreTest(
    int blockSize) {
  this->work = new MultipleTraitLinearRegressionScoreTestInternal;
  this->blockSize = blockSize;
  this->resultLength = 0;
  this->groupSize = -1;
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
  w.N = pheno.rows;
  w.T = pheno.cols;
  w.C = cov.cols;
  w.M = -1;

  w.Y.resize(tests.size());
  w.Z.resize(tests.size());
  w.ZZinv.resize(tests.size());
  w.hasCovariate.resize(tests.size());
  w.missingIndex.resize(tests.size());
  w.Uyz.resize(tests.size());
  w.Ugz.resize(tests.size());
  w.Uyg.resize(tests.size());
  w.sigma2.resize(tests.size());
  w.nTest = tests.size();
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
  std::vector<std::vector<std::string> > allCovName;
  // arrange Y, Z according to missing pattern for each trait
  for (int i = 0; i < w.nTest; ++i) {
    phenoName = tests.getPhenotype(i);
    phenoCol.clear();
    phenoCol.push_back(phenoDict[phenoName[0]]);
    covName = tests.getCovariate(i);
    allCovName.push_back(covName);
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

    // create missing indicators for Y[i] (and Z[i] if covariates used)
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
    if (w.Y[i].rows() == 0) {
      fprintf(stderr, "Due to missingness, there is no sample to test!\n");
      return -1;
    }

    // center and scale Y, Z
    scale(&w.Y[i]);
    scale(&w.Z[i]);

    // calcualte Uzy, inv(Z'Z)
    if (w.hasCovariate[i]) {
      w.ZZinv[i].noalias() =
          (w.Z[i].transpose() * w.Z[i])
              .ldlt()
              .solve(EMat::Identity(w.Z[i].cols(), w.Z[i].cols()));
      w.Uyz[i].noalias() = w.Z[i].transpose() * w.Y[i];
      w.sigma2[i] = (w.Y[i].transpose() * w.Y[i] -
                     w.Uyz[i].transpose() * w.ZZinv[i] * w.Uyz[i])(0, 0) /
                    w.Y[i].rows();
    } else {
      w.sigma2[i] = w.Y[i].col(0).squaredNorm() / w.Y[i].rows();
    }
  }  // end for i

  // Make groups based on model covariats and missing patterns of (Y, Z)
  // Detail:
  // For test: 1, 2, 3, ..., nTest, a possible grouping is:
  // (1, 3), (2), (4, 5) ...
  // =>
  // test_1 => group 0, offset 0
  // test_2 => group 1, offset 0
  // test_3 => group 0, offset 1
  //
  // For each test, we will use
  // [covar_name_1, covar_name_2, ...., missing_pattern],
  // as the dict key to distingish groups, and the dict value is the index of
  // the test
  std::map<std::vector<std::string>, int> groupDict;
  groupSize = 0;
  for (int i = 0; i < w.nTest; ++i) {
    std::vector<std::string> key = allCovName[i];
    key.push_back(toString(w.missingIndex[i]));
    if (0 == groupDict.count(key)) {
      groupDict[key] = groupSize;
      group.resize(groupSize + 1);
      group[groupSize].push_back(i);
      groupSize++;
    } else {
      group[groupDict[key]].push_back(i);
    }
  }
  // fprintf(stderr, "total %d missingness group\n", groupSize);

  w.G.resize(groupSize);
  w.groupedY.resize(groupSize);
  w.groupedZ.resize(groupSize);
  w.groupedUyz.resize(groupSize);
  w.groupedZZinv.resize(groupSize);
  w.groupedL.resize(groupSize);
  w.ustat.resize(groupSize);
  w.vstat.resize(groupSize);
  w.groupedHasCovariate.resize(groupSize);
  for (int i = 0; i < groupSize; ++i) {
    const int nc = group[i].size();
    const int nr = w.Y[group[i][0]].rows();
    w.groupedY[i].resize(nr, nc);
    for (int j = 0; j < nc; ++j) {
      w.groupedY[i].col(j) = w.Y[group[i][j]];
    }
    // initialize G
    w.G[i].resize(nr, blockSize);
    w.ustat[i].resize(blockSize, nc);
    w.vstat[i].resize(blockSize, 1);
    w.groupedZ[i] = w.Z[group[i][0]];
    w.groupedHasCovariate[i] = w.hasCovariate[group[i][0]];
    if (w.groupedHasCovariate[i]) {
      w.groupedUyz[i] = w.groupedZ[i].transpose() * w.groupedY[i];
    }
    w.groupedZZinv[i] = w.ZZinv[group[i][0]];
    Eigen::LLT<Eigen::MatrixXf> lltOfA(w.groupedZZinv[i]);
    // L * L' = A
    w.groupedL[i] = lltOfA.matrixL();

    // fprintf(stderr, "i = %d, group has covar = %s\n", i,
    //         w.groupedHasCovariate[i] ? "true" : "false");
  }
  // clean up memory
  w.Y.clear();
  w.Z.clear();
  w.Uyz.clear();
  w.hasCovariate.clear();
  w.ZZinv.clear();

  return true;
}

bool MultipleTraitLinearRegressionScoreTest::AddGenotype(const Matrix& g) {
  MultipleTraitLinearRegressionScoreTestInternal& w = *this->work;
  assert(resultLength < blockSize);
  for (int i = 0; i < groupSize; ++i) {
    // Assign g to G[, resultLength] accoring to missing pattern
    const std::vector<bool>& missingIndex =
        w.missingIndex[group[i][0]];  // 0: within the same group, missingIndex
                                      // are the same
    const int n = missingIndex.size();
    EMat& G = w.G[i];

    int idx = 0;
    for (int j = 0; j < n; ++j) {
      if (!missingIndex[j]) {
        G(idx, resultLength) = g[j][0];
        ++idx;
      }
    }
    assert(idx == G.rows());
  }
  resultLength++;
  return true;
}

bool MultipleTraitLinearRegressionScoreTest::TestCovariateBlock() {
  MultipleTraitLinearRegressionScoreTestInternal& w = *this->work;
  for (int i = 0; i < groupSize; ++i) {
    // delcare const variables
    EMat& G = w.G[i];
    EMat& Ugz = w.Ugz[i];
    EMat& Uyg = w.Uyg[i];
    const EMat& Z = w.groupedZ[i];
    const EMat& Y = w.groupedY[i];
    const EMat& Uyz = w.groupedUyz[i];
    const bool& hasCovariate = w.groupedHasCovariate[i];
    const EMat& ZZinv = w.groupedZZinv[i];
    const EMat& L = w.groupedL[i];

    // center and scale g
    scale(&G);

    // calculate Ugz, Uyg
    if (hasCovariate) {
      Ugz.noalias() = Z.transpose() * G;  // C by blockSize
    }
    Uyg.noalias() = G.transpose() * Y;  // blockSize by T

    // calculate Ustat, Vstat
    if (hasCovariate) {
      w.ustat[i].noalias() =
          (Uyg - Ugz.transpose() * ZZinv * Uyz);  // blockSize by T
      w.vstat[i].noalias() =
          ((G.array().square()).matrix().colwise().sum() -
           (L.transpose() * Ugz).array().square().matrix().colwise().sum())
              .transpose();  // blockSize by 1
    } else {                 // no covariate
      w.ustat[i].noalias() = Uyg;
      w.vstat[i].noalias() =
          G.array().square().matrix().colwise().sum().transpose();  // blockSize
                                                                    // by 1
    }
    // defer this, as this cannot be grouped
    // w.vstat[i] *= w.sigma2[group[i][0]];
  }
  // assign u, v; calculat p-values
  for (int j = 0; j < blockSize; ++j) {
    for (int i = 0; i < groupSize; ++i) {
      for (size_t k = 0; k != group[i].size(); ++k) {
        const int idx = group[i][k];
        const double u = w.ustat[i](j, k);
        const double v = w.vstat[i](j, 0) * w.sigma2[idx];
        this->ustat[j][idx] = u;
        this->vstat[j][idx] = v;

        if (v == 0.) {
          this->pvalue[j][idx] = NAN;
        } else {
          double stat = u * u / v;
          this->pvalue[j][idx] = gsl_cdf_chisq_Q(stat, 1.0);
        }
      }
    }
  }  // end for i
  return true;
}
