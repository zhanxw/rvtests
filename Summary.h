#ifndef _SUMMARY_H_
#define _SUMMARY_H_

#include "CommonFunction.h"

class Summary{
public:
Summary(): min(0), q1(0), median(0), q3(0), max(0), mean(0), sd(0), n(0){};
  void add(const std::vector<double>& v) {
    n = v.size();

    std::vector<double> t = v;
    std::sort(t.begin(), t.end());

    min = t[0];
    q1 = t[ (int) n * 0.25];
    median = t[ (int) n * 0.5];
    q3 = t[ (int) n * 0.75];
    max = t[ n - 1];

    mean = calculateMean(v);
    sd = calculateSD(v);
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
  void recordRawPhenotype(const std::vector<double>& phenoInOrder){
    rawPheno.add(phenoInOrder);
  };
  void recordTransformedPhenotype(const std::vector<double>& phenoInOrder, bool inverseNormalized){
    this->inverseNormal = inverseNormalized;
    transformedPheno.add(phenoInOrder);
  };
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
    int nr = m.rows;
    int nc = m.cols;
    this->covLabel.clear();
    this->cov.clear();
    for (int i = 0 ; i < nc; ++i) {
      this->covLabel.push_back( m.GetColumnLabel(i));
      this->recordCovariateColumn(m, i );
    }
  };
  void outputHeader(FILE* fp) {
    // write summaries
    int nSample = rawPheno.n;
    fprintf(fp, "##Samples=%d\n", nSample);
    fprintf(fp, "##AnalyzedSamples=%d\n", nSample);
    fprintf(fp, "##Families=%d\n", nSample);
    fprintf(fp, "##AnalyzedFamilies=%d\n", nSample);
    fprintf(fp, "##Founders=%d\n", nSample);
    fprintf(fp, "##AnalyzedFounders=%d\n", nSample);
    fprintf(fp, "##InverseNormal=%s\n", inverseNormal ? "ON" : "OFF");
    // write summaries
    fprintf(fp, "##TraitSummary\tmin\t25th\tmedian\t75th\tmax\tmean\tvariance\n");
    fprintf(fp, "##Triat\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
            rawPheno.min,
            rawPheno.q1,
            rawPheno.median,
            rawPheno.q3,
            rawPheno.max,
            rawPheno.mean,
            rawPheno.sd);
    fprintf(fp, "##AnalyzedTrait\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
            transformedPheno.min,
            transformedPheno.q1,
            transformedPheno.median,
            transformedPheno.q3,
            transformedPheno.max,
            transformedPheno.mean,
            transformedPheno.sd);

    // write covariate
    fprintf(fp, "##Covariates=");
    for (size_t i = 0; i < cov.size(); ++i ) {
      fputs(covLabel[i].c_str(), fp);
      if (i)
        fputc(',', fp);
    }
    fputc('\n', fp);
    
    fprintf(fp, "##CovariateSummary\tmin\t25th\tmedian\t75th\tmax\tmean\tvariance\n");
    for (size_t i = 0; i < cov.size(); ++i ) {
      fprintf(fp, "##%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
              covLabel[i].c_str(),
              cov[i].min,
              cov[i].q1,
              cov[i].median,
              cov[i].q3,
              cov[i].max,
              cov[i].mean,
              cov[i].sd);
    }
  }
private:
  Summary rawPheno;
  Summary transformedPheno;
  bool inverseNormal;
  std::vector<std::string> covLabel;
  std::vector<Summary> cov;
};

#endif /* _SUMMARY_H_ */
