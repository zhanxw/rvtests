#ifndef _PERMUTATION_H_
#define _PERMUTATION_H_

#include "Result.h"

#if 0
class AdaptivePermutationCheck{
 public:
  void addStat(double d) {
    this->stats.push_back(d);
  }
  /// compare to @param s, see we can early stop
  /// Algorithm:
  /// use existsing permuted stats to obtain a range of 95% confidence interval
  bool needEarlyStop(double s) {
    double m, sd;
    getMeanAndSD(&m, &sd);
    double lb = m - 1.96 * sd;
    double ub = m + 1.96 * sd;
    if ( s < lb ) return true;
    return false;
  };
 private:
  double getMean() {
    double s = 0;
    for (int i = 0; i < stats.size(); i++) {
      s+= this->stats[i];
    }
  };
  void getMeanAndSD(double* mean, double* sd) {
    if (stats.size() <= 1) {
      *mean = 0.0;
      *sd = 0.0;
      return;
    }
    *mean = getMean();
    double s = 0.0;
    for (int i = 0; i < stats.size(); ++i) {
      double tmp = (stats[i] - (*mean));
      s += (tmp*tmp);
    }
    *sd = sqrt( s / (stats.size() - 1));
  };
  std::vector<double> stats;
};
#endif

class Permutation{
 public:
  Permutation():numPerm(10000), alpha(0.05) {};
  Permutation(int nPerm, double alpha):numPerm(nPerm), alpha(alpha),
                                       obs(-1.), actualPerm(-1), threshold(-1.), numX(-1), numEqual(-1){
    result.addHeader("NumPerm");
    result.addHeader("ActualPerm");
    result.addHeader("Stat");
    result.addHeader("NumGreater");
    result.addHeader("NumEqual");
    result.addHeader("PermPvalue");
  };
  /**
   * @param observation: observed statistics
   */
  void init(double observation) {
    this->obs = observation;
    this->actualPerm = 0;
    this->threshold = 1.0 * this->numPerm * this->alpha * 2;
    this->numX = 0;
    this->numEqual = 0;

  };
  /**
   * @return true if need more permutations
   */
  bool next() {
    if (this->actualPerm >= this->numPerm) return false;
    if (numX + numEqual > threshold){
      return false;
    }
    return true;
  };
  void add(double s) {
    this->actualPerm++;
    if ( s > this->obs) {
      numX ++;
    }
    if ( s == this->obs) {
      numEqual ++;
    }
  };
  double getPvalue() const {
    if (this->actualPerm == 0) return 1.0;
    return  1.0 * (this->numX + 0.5 * this->numEqual) / this->actualPerm;
  };
  void reset() {
    obs = 0.0;
    actualPerm = 0;
    threshold = 0;
    numX = 0;
    numEqual = 0;
  };
  void writeHeader(FileWriter* fp){
    /* fprintf(fp, "%s\t%s\t%s\t%s\t%s\t%s", */
    /*         "NumPerm", "ActualPerm", "Stat", "NumGreater", "NumEqual", "PermPvalue"); */
    result.writeHeader(fp);
  }
  void writeHeaderTab(FileWriter* fp){
    /* fprintf(fp, "%s\t%s\t%s\t%s\t%s\t%s", */
    /*         "NumPerm", "ActualPerm", "Stat", "NumGreater", "NumEqual", "PermPvalue"); */
    result.writeHeaderTab(fp);
  }
  void writeHeaderLine(FileWriter* fp){
    /* fprintf(fp, "%s\t%s\t%s\t%s\t%s\t%s", */
    /*         "NumPerm", "ActualPerm", "Stat", "NumGreater", "NumEqual", "PermPvalue"); */
    result.writeHeaderLine(fp);
  }
  void updateValue() {
    result.updateValue("NumPerm", this->numPerm);
    result.updateValue("ActualPerm", this->actualPerm);
    result.updateValue("Stat", this->obs);
    result.updateValue("NumGreater", this->numX);
    result.updateValue("NumEqual", this->numEqual);
    result.updateValue("PermPvalue", this->getPvalue());
  }

  void writeOutput(FileWriter* fp) {
    /* fprintf(fp, "%d\t%d\t%g\t%d\t%d\t%g", */
    /*         this->numPerm, */
    /*         this->actualPerm, */
    /*         this->obs, */
    /*         this->numX, */
    /*         this->numEqual, */
    /*         this->getPvalue()); */
    updateValue();
    result.writeValue(fp);
  }
  void writeOutputLine(FileWriter* fp) {
    updateValue();
    result.writeValueLine(fp);
  }


 private:
  int numPerm;
  double alpha;
  double obs;
  int actualPerm;
  int threshold;
  int numX;
  int numEqual;
  Result result;
}; // class Permutation



#endif /* _PERMUTATION_H_ */
