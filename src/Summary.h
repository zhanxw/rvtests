#ifndef _SUMMARY_H_
#define _SUMMARY_H_

#include <algorithm>
#include <vector>
#include "base/CommonFunction.h"  // calculateMean, calculateSampleSD
#include "base/IO.h"
#include "base/MathVector.h"
#include "base/SimpleMatrix.h"

inline void assign(const std::vector<double>& in, Vector* out) {
  if (!out) return;
  int n = in.size();
  out->Dimension(n);
  for (int i = 0; i < n; ++i) {
    (*out)[i] = in[i];
  }
}

extern const char* VERSION;

class Summary {
 public:
  Summary() : min(0), q1(0), median(0), q3(0), max(0), mean(0), sd(0), n(0) {}
  void add(const std::vector<double>& v) {
    n = v.size();
    if (n == 0) return;

    std::vector<double> t = v;
    std::sort(t.begin(), t.end());

    min = t[0];
    q1 = t[(int)n * 0.25];
    median = t[(int)n * 0.5];
    q3 = t[(int)n * 0.75];
    max = t[n - 1];

    mean = calculateMean(v);
    sd = calculateSampleSD(v);
  }

 public:
  double min;
  double q1;
  double median;
  double q3;
  double max;
  double mean;
  double sd;
  int n;
};

/**
 * A class to summarize phenotype and genotypes
 */
class SummaryHeader {
 public:
  SummaryHeader() : inverseNormalized(false) {}
  void recordPhenotype(const char* label, const std::vector<double>& pheno) {
    this->phenoLabel.push_back(label);
    Summary s;
    s.add(pheno);
    this->pheno.push_back(s);
  }
  void setInverseNormalize(bool b) { this->inverseNormalized = b; }

  void recordCovariateColumn(const SimpleMatrix& m, int col) {
    Summary s;
    s.add(m.extractCol(col));
    this->cov.push_back(s);
  }
  // void recordCovariateColumn(Matrix& m, int col) {
  //   int nr = m.rows;

  //   std::vector<double> v(nr);
  //   for (int i = 0; i < m.rows; ++i) {
  //     v[i] = m[i][col];
  //   }

  //   Summary s;
  //   s.add(v);
  //   this->cov.push_back(s);
  // }

  void recordCovariate(const SimpleMatrix& m) {
    int nc = m.ncol();
    this->covLabel = m.getColName();
    this->cov.clear();
    for (int i = 0; i < nc; ++i) {
      this->recordCovariateColumn(m, i);
    }
  }
  void recordEstimation(const SimpleMatrix& est) {
    this->fittedResidualModel = est;
  }
  // void recordCovariate(Matrix& m) {
  //   // int nr = m.rows;
  //   int nc = m.cols;
  //   this->covLabel.clear();
  //   this->cov.clear();
  //   for (int i = 0; i < nc; ++i) {
  //     this->covLabel.push_back(m.GetColumnLabel(i));
  //     this->recordCovariateColumn(m, i);
  //   }
  // }
  void outputHeader(FileWriter* fp) {
    // write summaries
    int nSample = pheno.size() ? pheno[0].n : 0;
    fp->printf("##ProgramName=Rvtests\n");
    fp->printf("##Version=%s\n", VERSION);
    fp->printf("##Samples=%d\n", nSample);
    fp->printf("##AnalyzedSamples=%d\n", nSample);
    fp->printf("##Families=%d\n", nSample);
    fp->printf("##AnalyzedFamilies=%d\n", nSample);
    fp->printf("##Founders=%d\n", nSample);
    fp->printf("##AnalyzedFounders=%d\n", nSample);
    fp->printf("##InverseNormal=%s\n", this->inverseNormalized ? "ON" : "OFF");

    // write summaries
    fp->write("##TraitSummary\tmin\t25th\tmedian\t75th\tmax\tmean\tvariance\n");
    for (size_t i = 0; i < pheno.size(); ++i) {
      fp->printf("##%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", phenoLabel[i].c_str(),
                 pheno[i].min, pheno[i].q1, pheno[i].median, pheno[i].q3,
                 pheno[i].max, pheno[i].mean, pheno[i].sd * pheno[i].sd);
    }

    // write covariate
    if (!cov.empty()) {
      fp->write("##Covariates=");
      for (size_t i = 0; i < cov.size(); ++i) {
        if (i) fp->write(',');
        fp->write(covLabel[i].c_str());
      }
      fp->write('\n');

      fp->write(
          "##CovariateSummary\tmin\t25th\tmedian\t75th\tmax\tmean\tvariance\n");
      for (size_t i = 0; i < cov.size(); ++i) {
        fp->printf("##%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", covLabel[i].c_str(),
                   cov[i].min, cov[i].q1, cov[i].median, cov[i].q3, cov[i].max,
                   cov[i].mean, cov[i].sd * cov[i].sd);
      }
    }

    // optionally, write pheno ~ cov estimates
    // NOTE: when --useResidualAsPhenotype option is effective
    const int nParam = fittedResidualModel.nrow();
    if (nParam) {
      fp->write("##ResidualModelEstimates\n");
      fp->write("## - Name\tBeta\tSD\n");
      for (int i = 0; i < nParam; ++i) {
        fp->write("## - ");
        fp->write(fittedResidualModel.getRowName()[i]);
        fp->write("\t");
        if (finite(fittedResidualModel[i][0])) {
          fp->printf("%g", fittedResidualModel[i][0]);
        } else {
          fp->write("NA");
        }
        fp->write("\t");
        if (finite(fittedResidualModel[i][1])) {
          fp->printf("%g", fittedResidualModel[i][1]);
        } else {
          fp->write("NA");
        }
        fp->write("\n");
      }
    }
  }

  const std::vector<std::string>& getCovLabel() const { return this->covLabel; }

 private:
  std::vector<std::string> phenoLabel;
  std::vector<Summary> pheno;
  Summary transformedPheno;

  bool inverseNormalized;

  std::vector<std::string> covLabel;
  std::vector<Summary> cov;
  SimpleMatrix fittedResidualModel;
};

#endif /* _SUMMARY_H_ */
