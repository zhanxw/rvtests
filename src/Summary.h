#ifndef _SUMMARY_H_
#define _SUMMARY_H_

#include "base/IO.h"
#include "base/CommonFunction.h"
#include "ModelUtil.h"

inline void assign(const std::vector<double>& in, Vector* out) {
  if (!out) return;
  int n = in.size();
  out->Dimension(n);
  for (int i = 0; i < n; ++i) {
    (*out)[i] = in[i];
  }
}

extern const char* VERSION;

class Summary{
 public:
  Summary(): min(0), q1(0), median(0), q3(0), max(0), mean(0), sd(0), n(0){};
  void add(const std::vector<double>& v) {
    n = v.size();
    if (n == 0 ) return;

    std::vector<double> t = v;
    std::sort(t.begin(), t.end());

    min = t[0];
    q1 = t[ (int) n * 0.25];
    median = t[ (int) n * 0.5];
    q3 = t[ (int) n * 0.75];
    max = t[ n - 1];

    mean = calculateMean(v);
    sd = calculateSampleSD(v);
  };
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
class SummaryHeader{
 public:
  SummaryHeader(): inverseNormalized(false), isBinaryPhenotype(false), fitOK(false) {};
  void recordPhenotype(const char* label, const std::vector<double>& pheno){
    this->phenoLabel.push_back(label);
    Summary s;
    s.add(pheno);
    this->pheno.push_back(s);
  }
  void setInverseNormalize(bool b) {
    this->inverseNormalized = b;
  }

  /* void recordRawPhenotype(const std::vector<double>& phenoInOrder){ */
  /*   rawPheno.add(phenoInOrder); */
  /* }; */
  /* void recordTransformedPhenotype(const std::vector<double>& phenoInOrder, bool inverseNormalized){ */
  /*   this->inverseNormal = inverseNormalized; */
  /*   transformedPheno.add(phenoInOrder); */
  /* }; */
  void recordCovariateColumn(Matrix& m, int col) {
    int nr = m.rows;

    std::vector < double > v(nr);
    for (int i = 0; i < m.rows; ++i) {
      v[i] = m[i][col];
    }

    Summary s;
    s.add(v);
    this->cov.push_back(s);
  }
  void recordCovariate(Matrix& m) {
    // int nr = m.rows;
    int nc = m.cols;
    this->covLabel.clear();
    this->cov.clear();
    for (int i = 0 ; i < nc; ++i) {
      this->covLabel.push_back( m.GetColumnLabel(i));
      this->recordCovariateColumn(m, i );
    }
  }
  void fitModel(const std::vector<double>& pheno, bool binaryPhenotype, Matrix& cov) {
    this->isBinaryPhenotype = binaryPhenotype;
    Vector p;
    assign(pheno, &p);
    Matrix c;
    copyCovariateAndIntercept(p.Length(), cov, &c);
    if (binaryPhenotype) {
      fitOK = this->logistic.Fit(c, p);
    } else {
      fitOK = this->linear.Fit(c, p);
    }
  }
  void outputHeader(FileWriter* fp) {
    // write summaries
    int nSample = pheno.size()? pheno[0].n: 0;
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
      fp->printf("##%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
                 phenoLabel[i].c_str(),
                 pheno[i].min,
                 pheno[i].q1,
                 pheno[i].median,
                 pheno[i].q3,
                 pheno[i].max,
                 pheno[i].mean,
                 pheno[i].sd * pheno[i].sd );
    }

    if (!cov.empty()) {

      // write covariate
      fp->write("##Covariates=");
      for (size_t i = 0; i < cov.size(); ++i ) {
        if (i)
          fp->write(',');
        fp->write(covLabel[i].c_str());
      }
      fp->write('\n');

      fp->write("##CovariateSummary\tmin\t25th\tmedian\t75th\tmax\tmean\tvariance\n");
      for (size_t i = 0; i < cov.size(); ++i ) {
        fp->printf("##%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
                   covLabel[i].c_str(),
                   cov[i].min,
                   cov[i].q1,
                   cov[i].median,
                   cov[i].q3,
                   cov[i].max,
                   cov[i].mean,
                   cov[i].sd * cov[i].sd);
      }
    }

    //write null model
    fp->write("##NullModelEstimates\n");
    if (!fitOK) {
      fp->write("## - WARNING: Null model fit failed\n");
    } else {
      if (!this->isBinaryPhenotype) {
        printNullModelEstimate(fp,
                               linear.GetCovEst(),
                               linear.GetCovB(),
                               linear.GetSigma2());
      } else {
        printNullModelEstimate(fp,
                               logistic.GetCovEst(),
                               logistic.GetCovB(),
                               1.0);
      }
    }
  }
  void printNullModelEstimate(FileWriter* fp,
                              const Vector& beta,
                              const Matrix& betaSd,
                              const double sigma) {
    if (beta.Length() != betaSd.rows) {
      fprintf(stderr, "Dimension does not match in class Summary.\n");
      /* fprintf(stderr, "Dim %d %d %d\n", beta.Length(), betaSd.rows, (int)covLabel.size()); */
      return;
    }
    // NOTE: beta dimension and number of covariate labels may not equals,
    //       we may regression out covaraite of the response variables.
    /*     beta.Length() != (int)covLabel.size() + 1)  */
    fp->printf("## - Name\tBeta\tSD\n");
    
    // intercept
    fp->printf("## - Intercept\t%g\t%g\n", beta[0], betaSd[0][0]);
    // other cov
    const int n = covLabel.size();
    for (int i = 0; i < n; ++i) {
      if  (i + 1 >= beta.Length() ) break;
      fp->printf("## - %s\t%g\t%g\n", covLabel[i].c_str(), beta[i+1], betaSd[i+1][i+1]);
    }
    // sigma
    fp->printf("## - Sigma\t%g\tNA\n", sigma);
  }
 private:
  std::vector<std::string> phenoLabel;
  std::vector<Summary> pheno;
  Summary transformedPheno;

  bool inverseNormalized;

  std::vector<std::string> covLabel;
  std::vector<Summary> cov;

  bool isBinaryPhenotype;
  LinearRegression linear;
  LogisticRegression logistic;
  bool fitOK;
};

#endif /* _SUMMARY_H_ */
