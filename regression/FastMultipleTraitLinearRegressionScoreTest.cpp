#include "FastMultipleTraitLinearRegressionScoreTest.h"

#include <map>
#include <string>
#include <vector>

#include "gsl/gsl_cdf.h"  // use gsl_cdf_chisq_Q

#include "Eigen/Cholesky"  // ldlt
#include "Eigen/Dense"

#include "base/Utils.h"  // tolower
#include "regression/Formula.h"

typedef Eigen::VectorXf EVec;
typedef Eigen::MatrixXf EMat;
typedef std::vector<EMat> EMatVec;
typedef Eigen::Map<Eigen::MatrixXf> EMap;

// Stores column indice for y and z in the concatenated matrix (y||z)
struct TestIndex {
  int y;
  std::vector<int> z;
};

// Store sufficient statistics and scaling factors for each test
struct TestSet {
  float scale_xy;
  EMat xz;        // [ resultLength x C ]
  EVec scale_xz;  // [ C ]
  float scale_xx;
  EMat zz_inv;      // scaled
  EVec zy;          // scaled
  EMat L;           // derived from scaled zz_inv
  float sigma2;     // var(y)
  EVec af;          // block by 1
  EVec ustat;       // block by 1
  EVec vstat;       // block by 1
  EVec correction;  // block by 1, when x is rare, need to further correct it
  EVec indice;      // index for non-missing samples for given Y and Z
};

class FastMultipleTraitLinearRegressionScoreTestInternal {
 public:
  std::vector<TestIndex> testIndex;
  std::vector<TestSet> testSet;

  EMat Y_Z;    // N by (uniqT + uniqC)
  EMat yz;     // uniqT by uniqC (Y' * Z)
  EMat G;      // N by block
  EVec af;     // block by nTest
  EMat G_Y_Z;  // block by (uniqT + uniqC)
};

/// Column names of @param m are stored in @param dict
void addColNameToDict(Matrix& m, std::map<std::string, int>* dict) {
  std::map<std::string, int>& d = *dict;
  for (int i = 0; i < m.cols; ++i) {
    const int n = d.size();
    d[m.GetColumnLabel(i)] = n;
  }
}

#if 0
void createMissingInd(const EMat& m, EMat* out) {
  EMat& ret = *out;
  const int M = m.rows();
  const int N = m.cols();
  ret.resize(M, N);
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      ret(i, j) = std::isnan(m(i, j)) ? 1 : 0;
    }
  }
}
#endif

void createObsInd(const EMat& m, EMat* out) {
  EMat& ret = *out;
  const int M = m.rows();
  const int N = m.cols();
  ret.resize(M, N);
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      ret(i, j) = std::isnan(m(i, j)) ? 0 : 1;
    }
  }
}

/**
 * Scale matrix @param m to have zero mean and set missing values to be zero
 */
void center(EMat* m) {
  EMat& ret = *m;
  const int N = ret.rows();
  const int C = ret.cols();

  float sum;
  int nonMiss;
  float avg;
  for (int j = 0; j < C; ++j) {
    sum = 0.;
    nonMiss = 0;
    for (int i = 0; i < N; ++i) {
      if (std::isnan(ret(i, j))) {
        continue;
      }
      nonMiss++;
      sum += ret(i, j);
    }
    if (nonMiss) {
      avg = sum / nonMiss;
    } else {
      avg = 0.0;
    }
    for (int i = 0; i < N; ++i) {
      if (std::isnan(ret(i, j))) {
        ret(i, j) = 0;
      } else {
        ret(i, j) -= avg;
      }
    }
  }
  // this won't work - does not handle missing values.
  // (*m).rowwise() -= (*m).colwise().sum() / (*m).rows();
}

void extract(const EMat& in, int idx, EMat* out) { (*out) = in.col(idx); }
EMat extract(const EMat& in, int idx) { return in.col(idx); }

void extract(const EMat& in, const std::vector<int>& idx, EMat* out) {
  (*out).resize(in.rows(), idx.size());
  for (size_t i = 0; i != idx.size(); ++i) {
    (*out).col(i) = in.col(idx[i]);
  }
}
EMat extract(const EMat& in, const std::vector<int>& idx) {
  EMat out;
  extract(in, idx, &out);
  return out;
}

void extract(const EMat& in, const std::vector<int>& rowIdx, const int colIdx,
             EMat* out) {
  EMat& ret = *out;
  const int N = rowIdx.size();
  const int M = 1;
  ret.resize(N, M);
  for (int i = 0; i < N; ++i) {
    ret(i, 0) = in(rowIdx[i], colIdx);
  }
}

EMat extract(const EMat& in, const std::vector<int>& rowIdx, const int colIdx) {
  EMat out;
  extract(in, rowIdx, colIdx, &out);
  return out;
}

void extract(const EMat& in, const std::vector<int>& rowIdx,
             const std::vector<int>& colIdx, EMat* out) {
  EMat& ret = *out;
  const int N = rowIdx.size();
  const int M = colIdx.size();
  ret.resize(N, M);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      ret(i, j) = in(rowIdx[i], colIdx[j]);
    }
  }
}

EMat extract(const EMat& in, const std::vector<int>& rowIdx,
             const std::vector<int>& colIdx) {
  EMat out;
  extract(in, rowIdx, colIdx, &out);
  return out;
}

/**
 * @param in vector
 * @param x scalar
 * @return a vector = (in[0] - x, in[1] - x, ...)
 */
template <typename T>
std::vector<T> vectorMinus(const std::vector<T>& in, T x) {
  std::vector<T> out(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    out[i] = in[i] - x;
  }
  return out;
}

#if 0
void adjustForRareAllele(const EMat& g, const EVec& indice, int resultLen,
                         const EMat& vstat, Matrix* outAf, EVec* correctScale) {
  const int N = g.rows();
  assert(N == indice.size());
  assert(correctScale);
  assert(outAf);

  Matrix& af = *outAf;
  int threshold = sqrt(2.0 * N);
  int nModel = indice.sum();
  (*correctScale).resize(resultLen);
  for (int i = 0; i < resultLen; ++i) {
    float sumGeno = (g.col(i).transpose().array() * indice.array()).sum();
    af[i] = 0.5 * sumGeno / nModel;
    if (sumGeno < threshold) {
      (*correctScale)(i) = 2.0 * af * (1.0 - 2.0 * af) * nModel / vstat(i);
    } else {
      (*correctScale)(i) = 1.0;
    }
  }
}
#endif

FastMultipleTraitLinearRegressionScoreTest::
    FastMultipleTraitLinearRegressionScoreTest(int blockSize) {
  this->work = new FastMultipleTraitLinearRegressionScoreTestInternal;
  this->blockSize = blockSize;
  this->resultLength = 0;
}

FastMultipleTraitLinearRegressionScoreTest::
    ~FastMultipleTraitLinearRegressionScoreTest() {
  if (this->work) {
    delete this->work;
    this->work = NULL;
  }
}

bool FastMultipleTraitLinearRegressionScoreTest::FitNullModel(
    Matrix& cov, Matrix& pheno, const FormulaVector& tests) {
  FastMultipleTraitLinearRegressionScoreTestInternal& w = *this->work;
  // set some values
  const int N = pheno.rows;
  const int nTest = tests.size();

  std::vector<TestIndex>& testIndex = w.testIndex;
  std::vector<TestSet>& testSet = w.testSet;
  testIndex.resize(nTest);
  testSet.resize(nTest);

  // create dict (key: phenotype/cov name, val: index)
  std::map<std::string, int> pheno2Idx;
  addColNameToDict(pheno, &pheno2Idx);
  std::map<std::string, int> cov2Idx;
  addColNameToDict(cov, &cov2Idx);

  std::map<std::string, int> phenoDict;
  std::map<std::string, int> covDict;
  for (int i = 0; i < nTest; ++i) {
    const std::vector<std::string>& phenoName = tests.getPhenotype(i);
    const std::vector<std::string>& covName = tests.getCovariate(i);

    if (phenoDict.count(phenoName[0])) {
      // do nothing
    } else {
      phenoDict.insert(std::make_pair(phenoName[0], phenoDict.size()));
    }
    testIndex[i].y = phenoDict[phenoName[0]];

    for (size_t j = 0; j != covName.size(); ++j) {
      if (covName[j] == "1" || tolower(covName[j]) == "intercept") {
        continue;
      }
      if (covDict.count(covName[j])) {
        // do nothing
      } else {
        covDict.insert(std::make_pair(covName[j], covDict.size()));
      }
      testIndex[i].z.push_back(covDict[covName[j]]);
    }
  }

  const int uniqT = phenoDict.size();
  const int uniqC = covDict.size();
  for (int i = 0; i < nTest; ++i) {
    for (size_t j = 0; j < testIndex[i].z.size(); ++j) {
      testIndex[i].z[j] += uniqT;
    }
  }

  // w.Y_Z = cbind(uniq_pheno, uniq_cov)
  w.Y_Z.resize(N, uniqT + uniqC);
  for (std::map<std::string, int>::const_iterator iter = phenoDict.begin();
       iter != phenoDict.end(); ++iter) {
    const int phenoCol = pheno2Idx[iter->first];
    const int yzCol = iter->second;
    for (int i = 0; i < N; ++i) {
      w.Y_Z(i, yzCol) = pheno[i][phenoCol];
    }
  }
  for (std::map<std::string, int>::const_iterator iter = covDict.begin();
       iter != covDict.end(); ++iter) {
    const int covCol = cov2Idx[iter->first];
    const int yzCol = iter->second + uniqT;
    for (int i = 0; i < N; ++i) {
      w.Y_Z(i, yzCol) = cov[i][covCol];
    }
  }

  // create missing index
  EMat ind_Y_Z;
  createObsInd(w.Y_Z, &ind_Y_Z);

  // 1. center the matrix
  // 2. set missing values to zero
  center(&w.Y_Z);

  // Need to create zy for each model
  EMat zy = w.Y_Z.rightCols(uniqC).transpose() * w.Y_Z.leftCols(uniqT);

  // allocate memory
  for (int i = 0; i < nTest; ++i) {
    const int C = testIndex[i].z.size();
    TestSet& ts = testSet[i];
    // ts.xy.resize(blockSize);
    ts.xz.resize(blockSize, C);
    ts.scale_xz.resize(C);
    ts.zz_inv.resize(C, C);
    ts.zy.resize(C);
    ts.L.resize(C, C);
    ts.af.resize(blockSize);
    ts.ustat.resize(blockSize);
    ts.vstat.resize(blockSize);
    ts.correction.resize(blockSize);
  }

  // pre-calculate sufficient statistics
  for (int i = 0; i < nTest; ++i) {
    TestSet& ts = testSet[i];
    EVec indY = extract(ind_Y_Z, testIndex[i].y);
    EMat indZ = extract(ind_Y_Z, testIndex[i].z);
    EVec indModel = indY.array() * indZ.rowwise().prod().array();
    EVec indZY = indZ.transpose() * indY;

    const float OBS_MODEL = indModel.col(0).sum();
    ts.scale_xy = (OBS_MODEL / indY.col(0).sum());
    ts.scale_xz = OBS_MODEL / indZ.colwise().sum().array();
    ts.scale_xx = OBS_MODEL / N;
    ts.indice = indModel;

    EMat Y = extract(w.Y_Z, testIndex[i].y);
    EMat Z = extract(w.Y_Z, testIndex[i].z);

    ts.sigma2 = (Y.col(0).array() * indY.col(0).array()).square().sum() *
                OBS_MODEL / indY.sum();
    // handle covariate
    if (testIndex[i].z.size()) {
      ts.zz_inv = ((Z.transpose() * Z).array() * OBS_MODEL /
                   (indZ.transpose() * indZ).array())
                      .matrix()
                      .ldlt()
                      .solve(EMat::Identity(Z.cols(), Z.cols()));
      Eigen::LLT<Eigen::MatrixXf> lltOfA(ts.zz_inv);
      // L * L' = A
      ts.L = lltOfA.matrixL();

      ts.zy = extract(zy, vectorMinus(testIndex[i].z, uniqT), testIndex[i].y)
                  .array() *
              OBS_MODEL / indZY.array();

      ts.sigma2 -= (ts.zy.transpose() * ts.zz_inv * ts.zy)(0, 0);
    }
    ts.sigma2 /= OBS_MODEL;
  }

  // allocate memory for results
  w.G.resize(N, blockSize);
  this->ustat.Dimension(blockSize, tests.size());
  this->vstat.Dimension(blockSize, tests.size());
  this->pvalue.Dimension(blockSize, tests.size());

  return true;
}

bool FastMultipleTraitLinearRegressionScoreTest::AddGenotype(const Matrix& g) {
  EMat& G = this->work->G;
  assert(resultLength < blockSize);
  const int N = g.rows;
  for (int j = 0; j < N; ++j) {
    G(j, resultLength) = g[j][0];
  }

  resultLength++;
  return true;
}

bool FastMultipleTraitLinearRegressionScoreTest::TestCovariateBlock() {
  FastMultipleTraitLinearRegressionScoreTestInternal& w = *this->work;
  EMat& g = w.G;  // N by resultLength
  const float thresholdAC = sqrt(2.0 * g.rows());
  const int nTest = w.testSet.size();

  EVec nmiss;
  for (int i = 0; i < nTest; ++i) {
    TestSet& ts = w.testSet[i];
    // // correct V stat if necessary
    nmiss = g.transpose() * ts.indice;  // sum(geno, is.na = F)
    ts.af = 0.5 * nmiss /
            ts.indice.sum();  // num of non-missing elements in this test
    for (int j = 0; j < resultLength; ++j) {
      const float af = ts.af(j);
      if (nmiss(j) < thresholdAC) {
        ts.correction(j) = 2.0 * af * (1.0 - 2.0 * af) * ts.indice.sum();
      } else {
        ts.correction(j) = -1.0;
      }
    }
  }

  // calculate G'Y || G'Z
  // assume g does not have missing values
  center(&g);
  w.G_Y_Z = g.transpose() * w.Y_Z;  // resultLength x (uniqT + uniqC)

  for (int i = 0; i < nTest; ++i) {
    TestSet& ts = w.testSet[i];
    TestIndex& testIndex = w.testIndex[i];
    ts.ustat.noalias() =
        (extract(w.G_Y_Z, testIndex.y).array() * ts.scale_xy).matrix();
    ts.vstat.noalias() = ((g.colwise().squaredNorm()).transpose() *
                          ts.scale_xx);  // blockSize by 1

    for (int j = 0; j < resultLength; ++j) {
      if (ts.correction(j) > 0) {
        ts.correction(j) /= ts.vstat(j);
      } else {
        ts.correction(j) = 1.0;
      }
    }

    if (testIndex.z.size()) {
      ts.xz = extract(w.G_Y_Z, testIndex.z) * ts.scale_xz.asDiagonal();
      ts.ustat.noalias() -= ts.xz * ts.zz_inv * ts.zy;
      // g: [ N x resultLength ]
      // v = g' g  - g' z (z'z)^(-1) z' g
      //   = g' g -  (g' z L)' (L' z' g)
      // v_ii = norm(g.col(i) )^2 - norm((L' z' g).col(i))^2
      //      = norm(g.col(i) )^2 - norm((g' z L).row(i))^2
      ts.vstat.noalias() -=
          (ts.xz * ts.L).rowwise().squaredNorm();  // blockSize by 1
    }
    ts.vstat *= ts.sigma2;

    // apply correction
    ts.vstat.array() *= ts.correction.array();
  }

  // assign u, v; calculat p-values
  for (int i = 0; i < resultLength; ++i) {
    for (int j = 0; j < nTest; ++j) {
      TestSet& ts = w.testSet[j];
      const float u = ts.ustat(i);
      const float v = ts.vstat(i);
      this->ustat[i][j] = u;
      this->vstat[i][j] = v;

      if (this->vstat[i][j] == 0.0) {
        this->pvalue[i][j] = NAN;
      } else {
        float stat = u * u / v;
        this->pvalue[i][j] = gsl_cdf_chisq_Q(stat, 1.0);
      }
    }
  }
  return true;
}
