#include "ModelFitter.h"

#include "Model.h"
#include "ModelParser.h"

#include "regression/GSLIntegration.h"

//////////////////////////////////////////////////////////////////////
// Implementation of various collpasing methods
#if 0
/**
 * @return allele frequency
 * NOTE: Madson-Browning definition of alleleFrequency is different
 *       double freq = 1.0 * (ac + 1) / (an + 1);
 */
double getMarkerFrequency(Matrix& in, int col) {
  int& numPeople = in.rows;
  double ac = 0;  // NOTE: here genotype may be imputed, thus not integer
  int an = 0;
  for (int p = 0; p < numPeople; p++) {
    if ((int)in[p][col] >= 0) {
      ac += in[p][col];
      an += 2;
    }
  }
  if (an == 0) return 0.0;
  double freq = ac / an;
  return freq;
}

void getMarkerFrequency(Matrix& in, std::vector<double>* freq) {
  freq->resize(in.cols);
  for (int i = 0; i < in.cols; ++i) {
    (*freq)[i] = getMarkerFrequency(in, i);
  }
}
#endif
double getMarkerFrequency(DataConsolidator* dc, int col) {
  return dc->getMarkerFrequency(col);
}

void getMarkerFrequency(DataConsolidator* dc, std::vector<double>* freq) {
  return dc->getMarkerFrequency(freq);
}

double getMarkerFrequencyFromControl(Matrix& in, Vector& pheno, int col) {
  int& numPeople = in.rows;
  double ac = 0;  // NOTE: here genotype may be imputed, thus not integer
  int an = 0;
  for (int p = 0; p < numPeople; p++) {
    if (pheno[p] == 1) continue;
    if (in[p][col] >= 0) {
      ac += in[p][col];
      an += 2;
    }
  }
  // Refer:
  // 1. Madsen BE, Browning SR. A Groupwise Association Test for Rare Mutations
  // Using a Weighted Sum Statistic. PLoS Genet. 2009;5(2):e1000384. Available
  // at: http://dx.doi.org/10.1371/journal.pgen.1000384 [Accessed November 24,
  // 2010].
  double freq = 1.0 * (ac + 1) / (an + 2);
  return freq;
}

/**
 * Collapsing and combine method (indicator of existence of alternative allele)
 * @param in : sample by marker matrix
 * @param out: sample by 1 matrix
 */
void cmcCollapse(DataConsolidator* dc, Matrix& in, Matrix* out) {
  assert(out);
  int numPeople = in.rows;
  int numMarker = in.cols;

  out->Dimension(numPeople, 1);
  out->Zero();
  for (int p = 0; p < numPeople; p++) {
    for (int m = 0; m < numMarker; m++) {
      int g = (int)(in[p][m]);
      if (g > 0) {
        (*out)[p][0] = 1.0;
        break;
      }
    }
  }
}

void cmcCollapse(DataConsolidator* dc, Matrix& in,
                 const std::vector<int>& index, Matrix* out, int outIndex) {
  assert(out);
  int numPeople = in.rows;
  assert(out->rows == numPeople);
  assert(out->cols > outIndex);

  for (int p = 0; p < numPeople; p++) {
    (*out)[p][outIndex] = 0.0;
    for (size_t m = 0; m < index.size(); m++) {
      int g = (int)(in[p][index[m]]);
      if (g > 0) {
        (*out)[p][outIndex] = 1.0;
        break;
      }
    };
  };
}

/**
 * Morris-Zeggini method (count rare variants).
 * @param in : sample by marker matrix
 * @param out: sample by 1 matrix
 */
void zegginiCollapse(DataConsolidator* dc, Matrix& in, Matrix* out) {
  assert(out);
  int numPeople = in.rows;
  int numMarker = in.cols;

  out->Dimension(numPeople, 1);
  out->Zero();
  for (int p = 0; p < numPeople; p++) {
    for (int m = 0; m < numMarker; m++) {
      int g = (int)(in[p][m]);
      if (g > 0) {  // genotype is non-reference
        (*out)[p][0] += 1.0;
      }
    }
  }
}

void zegginiCollapse(DataConsolidator* dc, Matrix& in,
                     const std::vector<int>& index, Matrix* out, int outIndex) {
  assert(out);
  int numPeople = in.rows;
  assert(out->rows == numPeople);
  assert(out->cols > outIndex);

  for (int p = 0; p < numPeople; p++) {
    (*out)[p][outIndex] = 0.0;
    for (size_t m = 0; m < index.size(); m++) {
      int g = (int)(in[p][index[m]]);
      if (g > 0) {
        (*out)[p][outIndex] += 1.0;
      }
    };
  };
}

/**
 * @param genotype : people by marker matrix
 * @param phenotype: binary trait (0 or 1)
 * @param out: collapsed genotype
 */
void madsonBrowningCollapse(DataConsolidator* dc, Matrix& genotype,
                            Vector& phenotype, Matrix* out) {
  assert(out);
  int& numPeople = genotype.rows;
  int numMarker = genotype.cols;

  out->Dimension(numPeople, 1);
  out->Zero();

  for (int m = 0; m < numMarker; m++) {
    // calculate weight
    double freq = getMarkerFrequencyFromControl(genotype, phenotype, m);
    if (freq <= 0.0 || freq >= 1.0) continue;  // avoid freq == 1.0
    double weight = 1.0 / sqrt(freq * (1.0 - freq) * genotype.rows);
    // fprintf(stderr, "freq = %f\n", freq);

    for (int p = 0; p < numPeople; p++) {
      (*out)[p][0] += genotype[p][m] * weight;
    }
  }
};

void fpCollapse(DataConsolidator* dc, Matrix& in, Matrix* out) {
  assert(out);
  int& numPeople = in.rows;
  int numMarker = in.cols;

  out->Dimension(numPeople, 1);
  out->Zero();

  for (int m = 0; m < numMarker; m++) {
    // calculate weight
    // double freq = getMarkerFrequency(in, m);
    double freq = dc->getMarkerFrequency(m);
    if (freq <= 0.0 || freq >= 1.0) continue;  // avoid freq == 1.0
    double weight = 1.0 / sqrt(freq * (1.0 - freq));
    // fprintf(stderr, "freq = %f\n", freq);

    for (int p = 0; p < numPeople; p++) {
      (*out)[p][0] += in[p][m] * weight;
    }
  }
}

void madsonBrowningCollapse(DataConsolidator* dc, Matrix* d, Matrix* out) {
  assert(out);
  Matrix& in = (*d);
  int& numPeople = in.rows;
  int numMarker = in.cols;

  out->Dimension(numPeople, 1);
  out->Zero();

  for (int m = 0; m < numMarker; m++) {
    // calculate weight
    // double freq = getMarkerFrequency(in, m);
    double freq = dc->getMarkerFrequency(m);
    if (freq <= 0.0 || freq >= 1.0) continue;  // avoid freq == 1.0
    double weight = 1.0 / sqrt(freq * (1.0 - freq));
    // fprintf(stderr, "freq = %f\n", freq);

    for (int p = 0; p < numPeople; p++) {
      (*out)[p][0] += in[p][m] * weight;
    }
  };
}

/**
 * Convert genotype back to reference allele count
 * e.g. genotype 2 means homAlt/homAlt, so it has reference allele count 0
 */
/**
 * Convert genotype back to reference allele count
 * e.g. genotype 2 means homAlt/homAlt, so it has reference allele count 0
 */
void convertToReferenceAlleleCount(Matrix* g) {
  Matrix& m = *g;
  for (int i = 0; i < m.rows; ++i) {
    for (int j = 0; j < m.cols; ++j) {
      m[i][j] = 2 - m[i][j];
    }
  }
}

void convertToReferenceAlleleCount(Matrix& in, Matrix* g) {
  Matrix& m = *g;
  m = in;
  convertToReferenceAlleleCount(&m);
}

/**
 * group genotype by its frequency
 * @param in: sample by marker genotype matrix
 * @param out: key: frequency value:0-based index for freq
 * e.g. freq = [0.1, 0.2, 0.1, 0.3]  =>
 *      *group = {0.1: [0, 2], 0.2: 1, 0.3 : 3}
 * NOTE: due to rounding errors, we only keep 6 digits
 */
void groupFrequency(const std::vector<double>& freq,
                    std::map<double, std::vector<int> >* group) {
  group->clear();
  for (size_t i = 0; i != freq.size(); ++i) {
    double f = ceil(1000000. * freq[i]) / 1000000;
    (*group)[f].push_back(i);
  }
}

#if 0

/**
 * Collapsing @param in (people by marker) to @param out (people by marker),
 * if @param freqIn is empty, then frequncy is calculated from @param in
 * or according to @param freqIn to rearrange columns of @param in.
 * Reordered frequency are stored in @param freqOut, in ascending order
 */
void rearrangeGenotypeByFrequency(Matrix& in, const std::vector<double>& freqIn,
                                  Matrix* out, std::vector<double>* freqOut) {
  std::map<double, std::vector<int> > freqGroup;
  std::map<double, std::vector<int> >::const_iterator freqGroupIter;
  if (freqIn.empty()) {
    getMarkerFrequency(in, freqOut);
    groupFrequency(*freqOut, &freqGroup);
  } else {
    groupFrequency(freqIn, &freqGroup);
  }

  Matrix& sortedGenotype = *out;
  sortedGenotype.Dimension(in.rows, freqGroup.size());
  sortedGenotype.Zero();
  freqOut->clear();
  int idx = 0;
  for (freqGroupIter = freqGroup.begin(); freqGroupIter != freqGroup.end();
       freqGroupIter++) {
    freqOut->push_back(freqGroupIter->first);
    const std::vector<int>& cols = freqGroupIter->second;
    for (size_t j = 0; j != cols.size(); ++j) {
      for (int i = 0; i < in.rows; ++i) {
        sortedGenotype[i][cols[j]] += in[i][cols[j]];
      }
    }
    ++idx;
  }
}
#endif

void makeVariableThreshodlGenotype(
    DataConsolidator* dc, Matrix& in, const std::vector<double>& freqIn,
    Matrix* out, std::vector<double>* freqOut,
    void (*collapseFunc)(DataConsolidator*, Matrix&, const std::vector<int>&,
                         Matrix*, int)) {
  assert((int)freqIn.size() == in.cols);
  assert(freqIn.size());
  assert(out);
  assert(freqOut);

  std::map<double, std::vector<int> > freqGroup;
  std::map<double, std::vector<int> >::const_iterator freqGroupIter;

  groupFrequency(freqIn, &freqGroup);

  Matrix& sortedGenotype = *out;
  sortedGenotype.Dimension(in.rows, freqGroup.size());
  sortedGenotype.Zero();
  freqOut->resize(freqGroup.size());
  int idx = 0;
  std::vector<int> cumCols;
  for (freqGroupIter = freqGroup.begin(); freqGroupIter != freqGroup.end();
       freqGroupIter++) {
    (*freqOut)[idx] = freqGroupIter->first;
    const std::vector<int>& cols = freqGroupIter->second;
    for (size_t i = 0; i != cols.size(); ++i) {
      cumCols.push_back(cols[i]);
    }
    (*collapseFunc)(dc, in, cumCols, out, idx);

#if 0
    printf("In:\n");
    print(in);
    printf("Out:\n");
    print(*out);
#endif
    ++idx;
  }
}

double fIntegrand(double x, void* param) {
  if (x > 500 || x < -500) return 0.0;
  const double alpha = *((double*)param);
  const double tmp = exp(alpha + x);
  const double k = 1.0 / sqrt(2.0 * 3.1415926535897);
  const double ret = tmp / (1. + tmp) / (1. + tmp) * k * exp(-x * x * 0.5);
  // fprintf(stderr, "alpha = %g\tx = %g\t ret = %g\n", alpha, x, ret);
  return ret;
}

int obtainB(double alpha, double* out) {
  Integration i;
  gsl_function F;
  F.function = &fIntegrand;
  F.params = &alpha;
  if (i.integrate(F)) {
    fprintf(stderr, "Calculation of b may be inaccurate.\n");
    return -1;
  }
  *out = i.getResult();
  return 0;
}

int obtainB(float alpha, float* out) {
  double ret;
  obtainB((double)alpha, &ret);
  *out = ret;
  return 0;
}

void makeVariableThreshodlGenotype(
    DataConsolidator* dc, Matrix& in, Matrix* out, std::vector<double>* freqOut,
    void (*collapseFunc)(DataConsolidator*, Matrix&, const std::vector<int>&,
                         Matrix*, int)) {
  std::vector<double> freqIn;
  // getMarkerFrequency(in, &freqIn);
  dc->getMarkerFrequency(&freqIn);

  makeVariableThreshodlGenotype(dc, in, freqIn, out, freqOut, collapseFunc);
}

#if 0
void SingleVariantScoreTest::calculateConstant(Matrix& phenotype) {
  int nCase = 0;
  int nCtrl = 0;
  for (int i = 0; i < phenotype.rows; ++i) {
    if (phenotype[i][0] == 1) {
      ++nCase;
    } else if (phenotype[i][0] == 0) {
      ++nCtrl;
    }
  }
  double alpha;
  if (nCtrl > 0) {
    alpha = log(1.0 * nCase / nCtrl);
  } else {
    alpha = 500.;
  }
  obtainB(alpha, &this->b);
  fprintf(stderr, "alpha = %g, b = %g\n", alpha, b);  
}
#endif

void MetaScoreTest::MetaFamBinary::calculateB() {
  obtainB(this->alpha, &this->b);
  // fprintf(stderr, "alpha = %g, b = %g\n", alpha, b);
  return;
}

//////////////////////////////////////////////////
// MetaCov
//////////////////////////////////////////////////
class MetaCovBase {
 public:
  MetaCovBase() : needToFitNullModel(true) {}
  virtual ~MetaCovBase() {}
  virtual int FitNullModel(Matrix& genotype, DataConsolidator* dc) = 0;
  // virtual int transformGenotype(Genotype* out, DataConsolidator* dc) = 0;
  // virtual int calculateXX(const Genotype& x1, const Genotype& x2,
  //                         float* covXX) = 0;
  // virtual int calculateXZ(const Genotype& x, std::vector<float>* covXZ) = 0;
  virtual int transformGenotype(FloatMatrixRef& out, DataConsolidator* dc) = 0;
  virtual int calculateXX(FloatMatrixRef& x1, FloatMatrixRef& x2,
                          float* covXX) = 0;
  virtual int calculateXZ(FloatMatrixRef& inGeno, FloatMatrixRef& outXZ) = 0;

  virtual int calculateZZ(Matrix* covZZ) = 0;
  bool needToFitNullModel;
  bool hemiRegion;

 protected:
  Matrix cov;
};

class MetaCovFamQtl : public MetaCovBase {
 public:
  MetaCovFamQtl() : metaCov(FastLMM::SCORE, FastLMM::MLE) {}
  int FitNullModel(Matrix& genotype, DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& covariate = dc->getCovariate();
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);

    bool fitOK;
    if (!hemiRegion) {
      fitOK = (0 ==
               metaCov.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(),
                                    *dc->getKinshipSForAuto()));
      this->U = dc->getKinshipUForAuto();
      this->S = dc->getKinshipSForAuto();
    } else {
      if (!dc->hasKinshipForX()) {
        fitOK = false;
        return -1;
      }
      fitOK = (0 ==
               metaCov.FitNullModel(cov, phenotype, *dc->getKinshipUForX(),
                                    *dc->getKinshipSForX()));
      this->U = dc->getKinshipUForX();
      this->S = dc->getKinshipSForX();
    }
    if (fitOK) {
      needToFitNullModel = false;
      return 0;
    }
    return -1;
  }
  int transformGenotype(FloatMatrixRef& geno, DataConsolidator* dc) {
    if (!hemiRegion) {
      metaCov.TransformCentered(geno, *dc->getKinshipUForAuto(),
                                *dc->getKinshipSForAuto());
    } else {
      metaCov.TransformCentered(geno, *dc->getKinshipUForX(),
                                *dc->getKinshipSForX());
    }
    return 0;
  }
  int calculateXX(FloatMatrixRef& g1, FloatMatrixRef& g2, float* covXX) {
    // const int nSample = g1.size();
    metaCov.GetCovXX(g1, g2, *U, *S, covXX);
    // (*covXX) /= nSample;
    return 0;
  }
  int calculateXZ(FloatMatrixRef& g, FloatMatrixRef& covXZ) {
    metaCov.GetCovXZ(g, *U, *S, covXZ);
    return 0;
  }
  int calculateZZ(Matrix* covZZ) {
    metaCov.GetCovZZ(covZZ);
    return 0;
  }

 private:
  const EigenMatrix* U;
  const EigenMatrix* S;
  FastLMM metaCov;
};  // class MetaCovFamQtl

class MetaCovUnrelatedQtl : public MetaCovBase {
  int FitNullModel(Matrix& genotype, DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& covariate = dc->getCovariate();
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);

    bool fitOK = linear.FitLinearModel(cov, phenotype);
    if (!fitOK) {
      return -1;
    }

    sigma2 = linear.GetSigma2();
    return 0;
  }
  int transformGenotype(FloatMatrixRef& out, DataConsolidator* dc) {
    //     if (out->empty()) {
    //       return 0.;
    //     }

    //     float sum = 0.0;
    //     const size_t n = out->size();
    // #pragma omp parallel for
    //     for (size_t i = 0; i < n; ++i) {
    //       sum += out->at(i);
    //     }
    //     const float avg = sum / n;
    // #pragma omp parallel for
    //     for (size_t i = 0; i < n; ++i) {
    //       (*out)[i] -= avg;
    //     }
    REF_TO_EIGEN(out, g);
    g = g.rowwise() - g.colwise().mean();
    return 0;
  }
  int calculateXX(FloatMatrixRef& x1, FloatMatrixRef& x2, float* covXX) {
    assert(x1.nrow_ > 0 && x2.nrow_ > 0);
    //     const int n = x1.size();
    //     if (n == 0) return 0.0;

    //     float sum_ij = 0.0;  // sum of genotype[,i]*genotype[,j]

    // #pragma omp parallel for reduction(+ : sum_ij)
    //     for (int c = 0; c < n; ++c) {  // iterator each people
    //       // fprintf(using namespace std;err, "weight[%d] = %g\n", (int)c,
    //       // weight[int(c)]);
    //       // sum_ij += w1[c] * w2[c] / this->weight[(int)c];
    //       sum_ij += x1[c] * x2[c];
    //     }
    //     *covXX = sum_ij / this->sigma2;  //  / n;

    REF_TO_EIGEN(x1, x1E);
    REF_TO_EIGEN(x2, x2E);
    *covXX = x1E.col(0).dot(x2E.col(0)) / this->sigma2;
    return 0;
  }
  int calculateXZ(FloatMatrixRef& x, FloatMatrixRef& covXZ) {
    //     const int nc = this->cov.cols;
    //     (*covXZ).resize(nc);

    //     float sum = 0.0;
    //     const int n = x.size();
    //     if (n == 0) return 0;

    //     for (int c = 0; c < nc; ++c) {
    //       sum = 0.0;
    // #pragma omp parallel for reduction(+ : sum)
    //       for (int i = 0; i < n; ++i) {
    //         sum += x[i] * cov[i][c];
    //       }
    //       (*covXZ)[c] = sum / this->sigma2;
    //     }

    REF_TO_EIGEN(x, xE);
    REF_TO_EIGEN(covXZ, covXZ_E);

    Eigen::MatrixXf covE;
    G_to_Eigen(cov, &covE);  // TODO: change matrix datat type to Eigen

    covXZ_E = xE.transpose() * covE / sigma2;
    return 0;
  }
  int calculateZZ(Matrix* covZZ) {
    Matrix c;
    c = cov;
    centerMatrix(&c);
    Matrix ct;
    ct.Transpose(c);
    covZZ->Product(ct, c);
    covZZ->Multiply(1.0 / this->sigma2);
    return 0;
  }

 private:
  LinearRegression linear;
  float sigma2;
};  // end MetaCovUnrelatedQtl

class MetaCovFamBinary : public MetaCovBase {
 public:
  MetaCovFamBinary() : metaCov(FastLMM::SCORE, FastLMM::MLE) {}
  int FitNullModel(Matrix& genotype, DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& covariate = dc->getCovariate();
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);

    // calculate alpha, b
    int nCase = 0;
    int nCtrl = 0;
    for (int i = 0; i < phenotype.rows; ++i) {
      if (phenotype[i][0] == 1) {
        ++nCase;
      } else if (phenotype[i][0] == 0) {
        ++nCtrl;
      }
    }
    if (nCtrl > 0) {
      alpha = log(1.0 * nCase / nCtrl);
    } else {
      alpha = 500.;
    }
    calculateB();

    bool fitOK;
    if (!hemiRegion) {
      fitOK = (0 ==
               metaCov.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(),
                                    *dc->getKinshipSForAuto()));
      this->U = dc->getKinshipUForAuto();
      this->S = dc->getKinshipSForAuto();
    } else {
      if (!dc->hasKinshipForX()) {
        fitOK = false;
        return -1;
      }
      fitOK = (0 ==
               metaCov.FitNullModel(cov, phenotype, *dc->getKinshipUForX(),
                                    *dc->getKinshipSForX()));
      this->U = dc->getKinshipUForX();
      this->S = dc->getKinshipSForX();
    }
    if (fitOK) {
      needToFitNullModel = false;
      return 0;
    }
    return -1;
  }
  int transformGenotype(FloatMatrixRef& out, DataConsolidator* dc) {
    if (!hemiRegion) {
      metaCov.TransformCentered(out, *dc->getKinshipUForAuto(),
                                *dc->getKinshipSForAuto());
    } else {
      metaCov.TransformCentered(out, *dc->getKinshipUForX(),
                                *dc->getKinshipSForX());
    }
    return 0;
  }

  int calculateXX(FloatMatrixRef& g1, FloatMatrixRef& g2, float* covXX) {
    metaCov.GetCovXX(g1, g2, *U, *S, covXX);
    (*covXX) *= b * b;
    return 0;
  }
  int calculateXZ(FloatMatrixRef& g, FloatMatrixRef& covXZ) {
    metaCov.GetCovXZ(g, *U, *S, covXZ);
    // const int n = covXZ->size();
    // for (int i = 0; i < n; ++i) {
    //   (*covXZ)[i] *= b * b;
    // }
    REF_TO_EIGEN(covXZ, covXZ_E);
    covXZ_E *= (b * b);
    return 0;
  }
  int calculateZZ(Matrix* covZZ) {
    metaCov.GetCovZZ(covZZ);
    covZZ->Multiply(b * b);
    return 0;
  }

 private:
  void calculateB() {
    obtainB(this->alpha, &this->b);
    return;
  }

 private:
  float alpha;
  float b;

 private:
  const EigenMatrix* U;
  const EigenMatrix* S;
  FastLMM metaCov;
};  // end MetaCovFamBinary

class MetaCovUnrelatedBinary : public MetaCovBase {
  int FitNullModel(Matrix& genotype, DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& covariate = dc->getCovariate();
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);

    bool fitOK;
    if (this->needToFitNullModel || dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      // copyPhenotype(phenotype, &this->pheno);
      fitOK = logistic.FitLogisticModel(cov, phenotype, 100);
      if (!fitOK) return -1;
      needToFitNullModel = false;

      Vector& y = logistic.GetPredicted();
      this->nSample = y.Length();
      // this->weight = logistic.GetVariance();
      Vector& w = logistic.GetVariance();
      G_to_Eigen(w, &weight);
    }

    return 0;
  }
  int transformGenotype(FloatMatrixRef& out, DataConsolidator* dc) {
    // do not need to transform genotypes
    return 0;
  }
  // covX1X2 = x1' * W * x2
  int calculateXX(FloatMatrixRef& x1, FloatMatrixRef& x2, float* covX1X2) {
    // float& covXX = *covX1X2;
    // covXX = 0.0;
    // for (int i = 0; i < nSample; ++i) {
    //   covXX += x1[i] * this->weight[i] * x2[i];
    // }

    REF_TO_EIGEN(x1, x1E);
    REF_TO_EIGEN(x2, x2E);
    *covX1X2 = (x1E.array() * weight.array() * x2E.array()).sum();

    return 0;
  }
  // covXZ = g' W Z where Z = (z1, z2, ... , zp)
  int calculateXZ(FloatMatrixRef& x, FloatMatrixRef& covXZ) {
    // const int nCov = cov.cols;
    // (*covXZ).resize(nCov);
    // for (int c = 0; c < nCov; ++c) {
    //   float s = 0.;
    //   for (int i = 0; i < nSample; ++i) {
    //     s += x[i] * weight[i] * cov[i][c];
    //   }
    //   (*covXZ)[c] = s;
    // }
    REF_TO_EIGEN(x, xE);
    REF_TO_EIGEN(covXZ, covXZ_E);
    Eigen::MatrixXf covE;
    G_to_Eigen(cov, &covE);
    covXZ_E = xE.transpose() * weight.asDiagonal() * covE;
    return 0;
  }
  int calculateZZ(Matrix* argCovZZ) {
    Matrix& covZZ = *argCovZZ;
    const int nCov = cov.cols;
    covZZ.Dimension(nCov, nCov);
    for (int i = 0; i < nCov; ++i) {
      for (int j = 0; j <= i; ++j) {
        covZZ[i][j] = 0.0;
        for (int k = 0; k < nSample; ++k) {
          covZZ[i][j] += cov[k][i] * weight(k) * cov[k][j];
        }
        if (j != i) {
          covZZ[j][i] = covZZ[i][j];
        }
      }
    }
    return 0;
  }

 private:
  int nSample;
  // Vector pheno;
  LogisticRegression logistic;
  Eigen::VectorXf weight;
};  // end class MetaCovUnrelatedBinary

class MetaCovFamQtlBolt : public MetaCovBase {
 public:
  MetaCovFamQtlBolt() { fprintf(stderr, "MetaCovFamQtlBolt model started\n"); }
  int FitNullModel(Matrix& genotype, DataConsolidator* dc) {
    const std::string& fn = dc->getBoltGenotypeFilePrefix();

    // fit null model
    bool fitOK = 0 == bolt_.FitNullModel(fn, &dc->getPhenotype());
    if (!fitOK) return -1;
    needToFitNullModel = false;
    return 0;
  }
  int transformGenotype(FloatMatrixRef& out, DataConsolidator* dc) {
    // no need to transform
    return 0;
  }
  int calculateXX(FloatMatrixRef& x1, FloatMatrixRef& x2, float* covXX) {
    bolt_.GetCovXX(x1, x2, covXX);
    return 0;
  }
  int calculateXZ(FloatMatrixRef& x, FloatMatrixRef& covXZ) { return 0; }
  int calculateZZ(Matrix* covZZ) { return 0; }

 private:
  BoltLMM bolt_;
};  // class MetaCovFamQtlBolt

MetaCovTest::MetaCovTest(int windowSize)
    : model(NULL),
      modelAuto(NULL),
      modelX(NULL),
      useBolt(false),
      fitOK(false),
      useFamilyModel(false),
      isHemiRegion(false) {
  this->modelName = "MetaCov";
  this->indexResult = true;
  this->numVariant = 0;
  this->nSample = -1;
  this->fout = NULL;
  this->windowSize = windowSize;
  result.addHeader("CHROM");
  result.addHeader("START_POS");
  result.addHeader("END_POS");
  result.addHeader("NUM_MARKER");
  result.addHeader("MARKER_POS");
  result.addHeader("COV");
}
MetaCovTest::~MetaCovTest() {
  while (queue.size() > 0) {
    printCovariance(fout, queue, isBinaryOutcome());
    genoPool.deallocate(queue.front().geno);
    genoCovPool.deallocate(queue.front().covXZ);
    queue.pop_front();
  }
  if (modelAuto) {
    delete modelAuto;
    modelAuto = NULL;
  }
  if (modelX) {
    delete modelX;
    modelX = NULL;
  }
}
int MetaCovTest::fitWithGivenGenotype(Matrix& genotype, DataConsolidator* dc) {
  // Matrix& phenotype = dc->getPhenotype();
  // Matrix& covariate = dc->getCovariate();
  Result& siteInfo = dc->getResult();
  this->isHemiRegion = dc->isHemiRegion(0);

  if (genotype.cols != 1) {
    fitOK = false;
    return -1;
  }
  if (genotype.rows == 0) {
    fitOK = false;
    return -1;
  }
  this->useFamilyModel = dc->hasKinship();
  if (nSample < 0) {  // uninitialized
    // calculate variance of y
    nSample = genotype.rows;
    nCovariate = dc->getCovariate().cols + 1;  // intercept
    genoPool.setChunkSize(nSample);
    genoCovPool.setChunkSize(nCovariate);
  }
  if (nSample != genotype.rows) {
    fprintf(stderr, "Sample size changed at [ %s:%s ]\n",
            siteInfo["CHROM"].c_str(), siteInfo["POS"].c_str());
    fitOK = false;
    return -1;
  }

  loci.pos.chrom = siteInfo["CHROM"];
  loci.pos.pos = atoi(siteInfo["POS"]);

  // assign loci.geno, and
  // check if this is a monomorphic site, if so, just skip it.
  // skip monomorphic sites
  if (isMonomorphicMarker(genotype, 0)) {
    fitOK = false;
    return -1;
  }

  // fit null model
  if (isHemiRegion) {
    model = modelX;
    if (!model) {
      model = modelX = createModel(this->useFamilyModel, isBinaryOutcome());
      model->hemiRegion = true;
    }
  } else {
    model = modelAuto;
    if (!model) {
      model = modelAuto = createModel(this->useFamilyModel, isBinaryOutcome());
      model->hemiRegion = false;
    }
  }
  if (!model) return -1;

  if (model->needToFitNullModel || dc->isPhenotypeUpdated() ||
      dc->isCovariateUpdated()) {
    // copyCovariateAndIntercept(genotype.rows, covariate, &cov);
    fitOK = (0 == model->FitNullModel(genotype, dc));
    if (!fitOK) return -1;
    model->needToFitNullModel = false;
  }

  // assign genotype to loci.geno
  // loci.geno.resize(nSample);
  // for (int i = 0; i < nSample; ++i) {
  //   loci.geno[i] = genotype[i][0];
  // }
  assignGenotype(genotype, loci.geno);
  loci.covXZ = genoCovPool.allocate();
  if (!useBolt) {
    // model->transformGenotype(&loci.geno, dc);
    // model->calculateXZ(loci.geno, &loci.covXZ);
    // const int numCovariate = dc->getCovariate().cols;
    FloatMatrixRef x(genoPool.chunk(loci.geno), nSample, 1);
    FloatMatrixRef xz(genoCovPool.chunk(loci.geno), 1, nCovariate);
    model->transformGenotype(x, dc);
    if (nCovariate) {
      model->calculateXZ(x, xz);
    }
    if (model->needToFitNullModel || dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
      model->calculateZZ(&this->covZZ);
      CholeskyInverseMatrix(this->covZZ, &this->covZZInv);
    }
  }
  fitOK = true;
  return 0;
}  // fitWithGivenGenotype

void MetaCovTest::assignGenotype(Matrix& genotype, Genotype& genoIdx) {
  genoIdx = genoPool.allocate();
  float* p = genoPool.chunk(genoIdx);
  for (int i = 0; i < nSample; ++i) {
    p[i] = genotype[i][0];
  }
}

int MetaCovTest::printCovariance(FileWriter* fp,
                                 const std::deque<Loci>& lociQueue,
                                 bool binaryOutcome) {
  std::deque<Loci>::const_iterator iter = lociQueue.begin();

  const size_t numMarker = lociQueue.size();
  position.resize(numMarker);
  this->covXX.resize(numMarker);
  float covX1X2;
  for (int idx = 0; iter != lociQueue.end(); ++iter, ++idx) {
    position[idx] = iter->pos.pos;
    FloatMatrixRef frontGeno(genoPool.chunk(lociQueue.front().geno), nSample,
                             1);
    FloatMatrixRef frontGenoCov(genoCovPool.chunk(lociQueue.front().covXZ), 1,
                                nCovariate);
    if (!useBolt) {
      FloatMatrixRef iterGeno(genoPool.chunk(iter->geno), nSample, 1);
      // model->calculateXX(lociQueue.front().geno, iter->geno, &covX1X2);
      model->calculateXX(frontGeno, iterGeno, &covX1X2);
      // this->covXX[idx] = computeScaledXX(covX1X2, lociQueue.front().covXZ,
      //                                    iter->covXZ, this->covZZInv);
      FloatMatrixRef iterGenoCov(genoCovPool.chunk(iter->geno), 1, nCovariate);
      this->covXX[idx] =
          computeScaledXX(covX1X2, frontGenoCov, iterGenoCov, this->covZZInv);

    } else {
      // model->calculateXX(lociQueue.front().geno, iter->geno,
      // &this->covXX[idx]);
      FloatMatrixRef iterGeno(genoPool.chunk(iter->geno), nSample, 1);
      model->calculateXX(frontGeno, iterGeno, &covX1X2);
      this->covXX[idx] = covX1X2;
    }
  }

  result.updateValue("CHROM", lociQueue.front().pos.chrom);
  result.updateValue("START_POS", lociQueue.front().pos.pos);
  result.updateValue("END_POS", lociQueue.back().pos.pos);
  result.updateValue("NUM_MARKER", (int)numMarker);

  static std::string s;
  s.resize(0);

  appendToString(position, &s);
  result.updateValue("MARKER_POS", s);

  // divide n is by convention, no particular meaning.
  // will divide n for covXX, covXZ and covZZ
  const float scale = 1.0 / nSample;
  s.clear();
  appendToString(this->covXX, scale, &s);
  if (outputGwama || binaryOutcome) {
    s += ':';
    FloatMatrixRef covXZMat(genoCovPool.chunk(lociQueue.front().covXZ),
                            nCovariate, 1);
    appendToString(covXZMat, scale, &s);

    s += ':';
    appendToString(this->covZZ, scale, &s);
  }
  result.updateValue("COV", s);
  result.writeValueLine(fp);
  return 0;
}  // printCovariance

MetaCovBase* MetaCovTest::createModel(bool familyModel, bool binaryOutcome) {
  MetaCovBase* ret = NULL;
  if (this->useBolt) {
    if (binaryOutcome) {
      fprintf(stderr, "BoltLMM does not support binary outcomes! Exit...\n");
      exit(1);
    }
    ret = new MetaCovFamQtlBolt;
    return ret;
  }

  if (familyModel && !binaryOutcome) {
    ret = new MetaCovFamQtl;
  }
  if (familyModel && binaryOutcome) {
    ret = new MetaCovFamBinary;
  }
  if (!familyModel && !binaryOutcome) {
    ret = new MetaCovUnrelatedQtl;
  }
  if (!familyModel && binaryOutcome) {
    ret = new MetaCovUnrelatedBinary;
  }
  return ret;
}

#if 0
void appendHeritability(FileWriter* fp, const FastLMM& model) {
  return;

  // TODO: handle empiricalkinship and pedigree kinship better
  // we estimate sigma2_g * K, but a formal defiinte is 2 * sigma2_g *K,
  // so multiply 0.5 to scale it.
  const double sigma2_g = model.GetSigmaG2() * 0.5;
  const double sigma2_e = model.GetSigmaE2();
  const double herit = (sigma2_g + sigma2_e == 0.) ? 0 : sigma2_g / (sigma2_g + sigma2_e);

  fp->printf("#Sigma2_g\t%g\n", sigma2_g);
  fp->printf("#Sigma2_e\t%g\n", sigma2_e);
  fp->printf("#Heritability\t%g\n", herit);
}

void appendHeritability(FileWriter* fp, const GrammarGamma& model) {
  return;

  // TODO: handle empiricalkinship and pedigree kinship better
  // we estimate sigma2_g * K, but a formal defiinte is 2 * sigma2_g *K,
  // so multiply 0.5 to scale it.
  const double sigma2_g = model.GetSigmaG2() * 0.5;
  const double sigma2_e = model.GetSigmaE2();
  const double herit = (sigma2_g + sigma2_e == 0.) ? 0 : sigma2_g / (sigma2_g + sigma2_e);

  fp->printf("#Sigma2_g\t%g\n", sigma2_g);
  fp->printf("#Sigma2_e\t%g\n", sigma2_e);
  fp->printf("#Heritability\t%g\n", herit);
}

#endif
