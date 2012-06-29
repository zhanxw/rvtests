#ifndef _MODELFITTER_H_
#define _MODELFITTER_H_

#include "libsrc/MathMatrix.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "regression/LogisticRegression.h"
#include "regression/LogisticRegressionScoreTest.h"
#include "regression/LogisticRegressionPermutationTest.h"
#include "regression/LinearRegression.h"
#include "regression/LinearRegressionScoreTest.h"
#include "regression/Skat.h"
#include "regression/Table2by2.h"
#include "regression/kbac_interface.h"

// various collapsing method
// they all take people by marker matrix
// and they won't take special care of missing genotypes
double getMarkerFrequency(Matrix& in, int col);
double getMarkerFrequencyFromControl(Matrix& in, Vector& pheno, int col);

void cmcCollapse(Matrix& in, Matrix* out);
void zegginiCollapse(Matrix& in, Matrix* out);
void fpCollapse(Matrix& in, Matrix* out);
void madsonBrowningCollapse(Matrix& genotype, Vector& phenotype, Matrix* out);

void rearrangeGenotypeByFrequency(Matrix& in, Matrix* out, std::vector<double>* freq);

void permute(Vector* v);
void centerVector(Vector* v);

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

class Permutation{
public:
Permutation():numPerm(10000), alpha(0.05) {};
Permutation(int nPerm, double alpha):numPerm(nPerm), alpha(alpha) {};  
  /**
   * @param observation: observed statistics
   */
  void init(double observation) {
    this->obs = observation;
    this->numPerm =  0;
    this->actualPerm = 0;
    this->threshold = this->numPerm * this->alpha * 2;
    this->numX = 0; 
    this->numEqual = 0;
  };
  /**
   * @return true if need more permutations
   */
  bool next() {
    if (this->actualPerm > this->numPerm) return false;
    if (numX + numX > threshold){
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
    return  0.5 * (this->numX + this->numEqual) / this->actualPerm;
  };
  void reset() {
  numPerm = 0;
  actualPerm = 0;
  threshold = 0;
  numX = 0;
  numEqual = 0;
  };
  void writeHeader(FILE* fp){
    fprintf(fp, "%s\t%s\t%s\t%s\t%s\t%s",
            "NumPerm", "ActualPerm", "STAT", "NumX", "NumEqual", "PermPvalue");
  }
  void writeOutput(FILE* fp) {
    fprintf(fp, "%d\t%d\t%g\t%d\t%d\t%g", 
            this->numPerm,
            this->actualPerm,
            this->obs,
            this->numX,
            this->numEqual,
            this->getPvalue());
  };
  
private:
  double obs;
  double alpha;
  int numPerm;
  int actualPerm;
  int threshold;
  int numX;
  int numEqual;
};

// take X, Y, Cov and fit model
// note, ModelFitter will use VCFData as READ-ONLY data structure,
// and collapsing results are stored internally.
class ModelFitter{
public:
  /* // fitting model */
  /* virtual int fit(Matrix& phenotype, Matrix& genotype) = 0; */
  // fitting model
  virtual int fit(Matrix& phenotype, Matrix& genotype, Matrix& covariate) = 0;
  // write result header
  virtual void writeHeader(FILE* fp, const char* prependString) = 0;
  // write model output
  virtual void writeOutput(FILE* fp, const char* prependString) = 0;

  ModelFitter(){
    this->modelName = "Unassigned_Model_Name";
    this->binaryOutcome = false; // default: using continuous outcome
  };
  const std::string& getModelName() const { return this->modelName; };
  // for particular class to call when fitting repeatedly
  // e.g. clear permutation counter
  // e.g. clear internal cache
  virtual void reset() {}; 
  // virtual void needFittingCovariate();
  bool isBinaryOutcome() const{
    return this->binaryOutcome;
  }
  void setBinaryOutcome() {
    this->binaryOutcome = true;
  };
  void setContinuousOutcome() {
    this->binaryOutcome = false;
  };
protected:
  std::string modelName;
  bool binaryOutcome;
}; // end ModelFitter

void copyPhenotype(Matrix& in, Vector* o){
  Vector& out = *o;
  out.Dimension(in.rows);
  for (int i = 0; i <in.rows; ++i){
    out[i] = in[i][0];
  }
};

/**
 * copy @param in to @param o, with first column being intercept
 */
void copyGenotypeWithIntercept(Matrix& in, Matrix* o){
  Matrix& out = *o;
  out.Dimension(in.rows, in.cols+1);
  for (int i = 0; i <in.rows; ++i){
    out[i][0] = 1.0;
    for (int j = 0; j < in.cols; ++j) {
      out[i][j + 1] = in[i][j];
    }
  }
  out.SetColumnLabel(0, "Intercept");
  for (int i = 0; i < in.cols; ++i)
    out.SetColumnLabel(i + 1, in.GetColumnLabel(i));
};
/**
 * copy vector of one, @param in and @param cov to @param o (essentially: o = cbind(1, in, cov))
 */
void copyGenotypeWithCovariateAndIntercept(Matrix& in, Matrix& cov, Matrix* o){
  Matrix& out = *o;
  out.Dimension(in.rows, 1+in.cols+cov.cols);

  for (int i = 0; i <in.rows; ++i){
    out[i][0] = 1.0;
  }
  out.SetColumnLabel(0, "Intercept");

  for (int j = 0; j < in.cols; ++j) {
    for (int i = 0; i <in.rows; ++i){
      out[i][1+j] = in[i][j];
    }
    out.SetColumnLabel(1+j, in.GetColumnLabel(j));
  }

  for (int j = 0; j < cov.cols; ++j) {
    for (int i = 0; i < cov.rows; ++i){
      out[i][1+j + in.cols] = cov[i][j];
    }
    out.SetColumnLabel(1+j + in.cols, cov.GetColumnLabel(j));
  }
};

/**
 * copy intercept and @param cov to @param o with its first column equaling to 1.0
 */
void copyCovariateAndIntercept(int n, Matrix& cov, Matrix* o){
  if (cov.cols == 0 ) {
    (*o).Dimension(n, 1);
    for (int i = 0; i < n; ++i) {
      (*o)[i][0] = 1.0;
    }
    return ;
  }
  (*o).Dimension(n, 1 + cov.cols);
  for (int i = 0; i < n; ++i) {
    (*o)[i][0] = 1.0;
    for (int j = 0; j <cov.cols; ++j){
      (*o)[i][j+1] = cov[i][j];
    }
  }
  return;
};

class SingleVariantWaldTest: public ModelFitter{
public:
  SingleVariantWaldTest(){
    this->modelName = "SingleWald";
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& cov) {
    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }
    copyPhenotype(phenotype, &this->Y);
    
    if (cov.cols) {
      copyGenotypeWithCovariateAndIntercept(genotype, cov, &this->X);
    } else {
      copyGenotypeWithIntercept(genotype, &this->X);
    }

    if (!isBinaryOutcome()) {
      fitOK = linear.FitLinearModel(this->X, this->Y);
    } else {
      fitOK = logistic.FitLogisticModel(this->X, this->Y, 100);
    }
    return (fitOK ? 0 : 1);
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fprintf(fp, "Test\tBeta\tSE\tPvalue\n");
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    // skip interecept (column 0)
    for (int i = 1; i < this->X.cols; ++i) {
      if (!fitOK) {
        fputs(prependString, fp);
        fprintf(fp, "%s\tNA\tNA\tNA\n", this->X.GetColumnLabel(i));
      } else {
        double beta, se, pval;
        if (!isBinaryOutcome()) {
          beta = linear.GetCovEst()[i];
          se = sqrt(linear.GetCovB()[i][i]);
          pval = linear.GetAsyPvalue()[i];
        } else {
          beta = logistic.GetCovEst()[i];
          se = sqrt(logistic.GetCovB()[i][i]);
          pval = logistic.GetAsyPvalue()[i];
        }
        fputs(prependString, fp);
        fprintf(fp, "%s\t%g\t%g\t%g\n", this->X.GetColumnLabel(i), beta, se, pval);
      }
    }
  };
private:
  Matrix X; // 1 + cov + geno
  Vector Y; // phenotype
  LinearRegression linear;
  LogisticRegression logistic;
  bool fitOK;
}; // SingleVariantWaldTest

class SingleVariantScoreTest: public ModelFitter{
public:
  SingleVariantScoreTest(){
    this->modelName = "SingleScore";
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& covariate) {
    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }
    nSample = genotype.rows;
    af = getMarkerFrequency(genotype, 0);

    Matrix cov;
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);

    copyPhenotype(phenotype, &this->pheno);
    if (!isBinaryOutcome()) {
      fitOK = linear.FitNullModel(cov, pheno);
      if (!fitOK) return -1;
      fitOK = linear.TestCovariate(cov, pheno, genotype);
    } else {
      fitOK = logistic.FitNullModel(cov, pheno, 100);
      if (!fitOK) return -1;
      fitOK = logistic.TestCovariate(cov, pheno, genotype);
    }
    return (fitOK ? 0 : -1);
  };

  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fprintf(fp, "AF\tStat\tDirection\tPvalue\n");
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (fitOK) {
      if (!isBinaryOutcome())
        fprintf(fp, "%g\t%g\t%c\t%g\n", af, linear.GetStat(), linear.GetU()[0][0] > 0 ? '+': '-', linear.GetPvalue());
      else
        fprintf(fp, "%g\t%g\t%c\t%g\n", af, logistic.GetStat(), logistic.GetU()[0][0] > 0 ? '+': '-', logistic.GetPvalue());
    }else
      fputs("NA\tNA\tNA\tNA\n", fp);
  };
private:
  double af;
  int nSample;
  Vector pheno;
  LinearRegressionScoreTest linear;
  LogisticRegressionScoreTest logistic;
  bool fitOK;
}; // SingleVariantScoreTest

class SingleVariantFisherExactTest: public ModelFitter{
public:
  SingleVariantFisherExactTest() {
    this->modelName = "FisherExact";
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fputs("Fisher.N00\t", fp);
    fputs("Fisher.N01\t", fp);
    fputs("Fisher.N10\t", fp);
    fputs("Fisher.N11\t", fp);
    fputs("CtrlAF\t", fp);
    fputs("CaseAF\t", fp);
    fputs("Fisher.PvalueTwoSide\t", fp);
    fputs("Fisher.PvalueLess\t", fp);
    fputs("Fisher.PvalueGreater\n", fp);
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& cov) {
    if (genotype.cols == 0 || !isBinaryOutcome()) {
      fitOK = false;
      return -1;
    };
    if (cov.cols ) {
      fitOK = false;
      return -1;
    };

    // fit model
    caseAC = 0;
    caseAN = 0;
    ctrlAC = 0;
    ctrlAN = 0;
    // step 1, fit two by two table
    int numPeople = genotype.rows;
    for (int i = 0; i < numPeople; i++) {
      int geno = genotype[i][0];
      int pheno = phenotype[i][0];
      if ( !(0 <= geno && geno <= 2) ) continue;
      if ( !(0 <= pheno && pheno <= 1)) continue;
      if (pheno == 1) {
        caseAC += geno;
        caseAN += 2;
      } else {
        ctrlAC += geno;
        ctrlAN += 2;
      }

      if (geno == 0)
        model.Increment(0, pheno);
      else
        model.Increment(1, pheno);
    }

    // step 2, calculate pvalue
    model.UpdateMarginSum();
    model.FullFastFisherExactTest();

    this->fitOK = true;
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (fitOK) {
      fprintf(fp, "%d\t", model.Get00());
      fprintf(fp, "%d\t", model.Get01());
      fprintf(fp, "%d\t", model.Get10());
      fprintf(fp, "%d\t", model.Get11());
      if (ctrlAN == 0) {
        fprintf(fp, "0\t");
      } else{
        fprintf(fp, "%g\t", 1.0 * ctrlAC / ctrlAN);
      }
      if (caseAN == 0) {
        fprintf(fp, "0\t");
      } else{
        fprintf(fp, "%g\t", 1.0 * caseAC / caseAN);
      }
      fprintf(fp, "%lf\t", model.getPExactTwoSided());
      fprintf(fp, "%lf\t", model.getPExactOneSidedLess());
      fprintf(fp, "%lf\n"  , model.getPExactOneSidedGreater());
    } else {
      fprintf(fp, "0\t0\t0\t0\tNA\tNA\tNA\n");
    }
  };
  void reset() {
    model.reset();
  };
private:
  Table2by2 model;
  int caseAC;
  int caseAN;
  int ctrlAC;
  int ctrlAN;

  bool fitOK;
}; // SingleVariantFisherExactTest

class CMCTest: public ModelFitter{
public:
  CMCTest() {
    this->modelName = "CMC";
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (isBinaryOutcome()) {
      fprintf(fp, "NumVariant\tCMC.Pvalue\n");
    } else {
      fprintf(fp, "NumVariant\tNonRefSite\tCMC.Pvalue\n");
    }
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& covariate) {
    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    Matrix cov;
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);
    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i< phenotype.rows; i++){
      pheno[i] = phenotype[i][0];
    }

    cmcCollapse(genotype, &collapsedGenotype);

    if (isBinaryOutcome()) {
      fitOK = logistic.FitNullModel(cov, pheno, 100);
      if (!fitOK) return -1;
      fitOK = logistic.TestCovariate(cov, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    } else {
      fitOK = linear.FitNullModel(cov, pheno);
      if (!fitOK) return -1;
      fitOK = linear.TestCovariate(cov, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    }
  };

  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (!fitOK) {
      fprintf(fp, "%d\tNA\n", this->numVariant);
    } else {
      if (isBinaryOutcome()) {
        fprintf(fp, "%d\t%f\n", this->numVariant, logistic.GetPvalue());
      } else {
        fprintf(fp, "%d\t%d\t%f\n", this->numVariant, this->totalNonRefSite(), linear.GetPvalue());
      }
    };
  };
private:
  int totalNonRefSite(){
    double s = 0.0;
    for (int i = 0; i < collapsedGenotype.rows; ++i) {
      s += collapsedGenotype[i][0];
    }
    return (int)(s);
  }

  Matrix collapsedGenotype;
  LogisticRegressionScoreTest logistic;
  LinearRegressionScoreTest linear;
  bool fitOK;
  int numVariant;
}; // CMCTest

class CMCFisherExactTest: public ModelFitter{
public:
  CMCFisherExactTest() {
    this->modelName = "CMCFisherExact";
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fputs("exactCMC.N00\t", fp);
    fputs("exactCMC.N01\t", fp);
    fputs("exactCMC.N10\t", fp);
    fputs("exactCMC.N11\t", fp);
    fputs("exactCMC.PvalueTwoSide\t", fp);
    fputs("exactCMC.PvalueLess\t", fp);
    fputs("exactCMC.PvalueGreater\n", fp);
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& cov) {
    if (genotype.cols == 0) {
      fitOK = false;
      return 1;
    };
    if (cov.cols ) {
      fitOK = false;
      return -1;
    };

    // collapsing
    cmcCollapse(genotype, &collapsedGenotype);

    // fit model
    // step 1, fit two by two table
    int numPeople = collapsedGenotype.rows;
    for (int i = 0; i < numPeople; i++) {
      int geno = collapsedGenotype[i][0];
      int pheno = phenotype[i][0];
      if ( !(0 <= geno && geno <= 1) ) continue;
      if ( !(0 <= pheno && pheno <= 1)) continue;
      model.Increment(geno, pheno);
    }

    // step 2, calculate pvalue
    model.UpdateMarginSum();
    model.FullFastFisherExactTest();

    this->fitOK = true;
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (fitOK) {
      fprintf(fp, "%d\t", model.Get00());
      fprintf(fp, "%d\t", model.Get01());
      fprintf(fp, "%d\t", model.Get10());
      fprintf(fp, "%d\t", model.Get11());

      fprintf(fp, "%lf\t", model.getPExactTwoSided());
      fprintf(fp, "%lf\t", model.getPExactOneSidedLess());
      fprintf(fp, "%lf\n"  , model.getPExactOneSidedGreater());
    } else {
      fprintf(fp, "0\t0\t0\t0\tNA\tNA\tNA\n");
    }
  };
  void reset() {
    model.reset();
  };
private:
  Matrix collapsedGenotype;
  Table2by2 model;
  bool fitOK;
}; // CMCFisherExactTest


class ZegginiTest: public ModelFitter{
public:
  ZegginiTest(){
    this->modelName = "ZegginiTest";
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fprintf(fp, "NumVariant\tZeggini.Pvalue\n");
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& covariate) {
    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    Matrix cov;
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);

    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i< phenotype.rows; i++){
      pheno[i] = phenotype[i][0];
    }

    zegginiCollapse(genotype, &collapsedGenotype);

    if (isBinaryOutcome()) {
      fitOK = logistic.FitNullModel(cov, pheno, 100);
      if (!fitOK) return -1;
      fitOK = logistic.TestCovariate(cov, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    } else {
      fitOK = linear.FitNullModel(cov, pheno);
      if (!fitOK) return -1;
      fitOK = linear.TestCovariate(cov, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    }
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (!fitOK) {
      fprintf(fp, "%d\tNA\n", this->numVariant);
    } else {
      if (isBinaryOutcome()) {
        fprintf(fp, "%d\t%f\n", this->numVariant, logistic.GetPvalue());
      } else {
        fprintf(fp, "%d\t%f\n", this->numVariant, linear.GetPvalue());
      }
    };
  };
private:
  Matrix collapsedGenotype;
  LogisticRegressionScoreTest logistic;
  LinearRegressionScoreTest linear;
  bool fitOK;
  int numVariant;
}; // ZegginiTest

class MadsonBrowningTest: public ModelFitter{
public:
MadsonBrowningTest(int nPerm, double alpha): perm(nPerm, alpha) {
    this->modelName = "MadsonBrowningTest";
  }
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& covariate) {
    if (!isBinaryOutcome()) {
      fitOK = false;
      return -1;
    }

    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    Matrix cov;
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);

    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i< phenotype.rows; i++){
      pheno[i] = phenotype[i][0];
    }

    madsonBrowningCollapse(genotype, pheno, &collapsedGenotype);


    fitOK = logistic.FitNullModel(cov, pheno, 100);
    if (!fitOK) return -1;
    fitOK = logistic.TestCovariate(cov, pheno, collapsedGenotype);
    if (!fitOK) return -1;

    // record observed stat
    this->perm.init(logistic.GetStat()); // a chi-dist
    
    while (this->perm.next()){
      int failed = 0;
      permute(&pheno);
      fitOK = logistic.TestCovariate(cov, pheno, collapsedGenotype);
      if (!fitOK) {
        if (failed < 10) {
          failed++;
          continue;
        }else {
          fitOK = false;
          return -1;
        }
      }
      // record new stats
      double pStat = logistic.GetStat();
      this->perm.add(pStat);
    } // end permutation
    return (fitOK ? 0 : -1);
  };
  void reset() {
    this->perm.reset();
  }
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (isBinaryOutcome()){
      fprintf(fp, "NumVariant\t");
      perm.writeHeader(fp);
      fprintf(fp, "\n");
    } else
      fprintf(fp, "NumVariant\tNA\n");
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fprintf(fp, "%d\tNA\n", this->numVariant);
    if (isBinaryOutcome()) {
      if (fitOK){
        perm.writeOutput(fp);
        fprintf(fp, "\n");
      } else {
        fprintf(fp, "NA\tNA\tNA\tNA\tNA\n");
      }
    } else {
      fprintf(fp, "%d\tNA\n", this->numVariant);
    }
  };
private:
  Matrix collapsedGenotype;
  LogisticRegressionScoreTest logistic;
  bool fitOK;
  int numVariant;
  Permutation perm;
  /* int nPerm; */
  /* int numX; */
  /* int actualPerm; */
  /* double stat; */
}; // MadsonBrowningTest

// Danyu Lin's method, using 1/sqrt(p(1-p)) as weight
// where p is estimated from all samples
class FpTest: public ModelFitter{
public:
  FpTest() {
    this->modelName = "FpTest";
  }
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fprintf(fp, "NumVariant\tFp.Pvalue\n");
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& covariate) {
    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    Matrix cov;
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);

    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i< phenotype.rows; i++){
      pheno[i] = phenotype[i][0];
    }

    fpCollapse(genotype, &collapsedGenotype);

    if (isBinaryOutcome()) {
      fitOK = logistic.FitNullModel(cov, pheno, 100);
      if (!fitOK) return -1;
      fitOK = logistic.TestCovariate(cov, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    } else {
      fitOK = linear.FitNullModel(cov, pheno);
      if (!fitOK) return -1;
      fitOK = linear.TestCovariate(cov, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    }
  };

  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (!fitOK) {
      fprintf(fp, "%d\tNA\n", this->numVariant);
    } else {
      if (isBinaryOutcome()) {
        fprintf(fp, "%d\t%f\n", this->numVariant, logistic.GetPvalue());
      } else {
        fprintf(fp, "%d\t%f\n", this->numVariant, linear.GetPvalue());
      }
    };
  };
private:
  Matrix collapsedGenotype;
  LogisticRegressionScoreTest logistic;
  LinearRegressionScoreTest linear;
  bool fitOK;
  int numVariant;
}; // FpTest

class RareCoverTest: public ModelFitter{
public:
RareCoverTest(int nPerm, double alpha): perm(nPerm, alpha) {
    this->modelName = "RareCover";
  }
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& covariate) {
    if (!isBinaryOutcome()) {
      fitOK = false;
      return -1;
    }
    if (covariate.cols == 0) {
      fitOK = false;
      return -1;
    };

    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    // use marker by people matrix for faster computation
    this->genotype.Transpose(genotype);
    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i< phenotype.rows; i++){
      pheno[i] = phenotype[i][0];
    }

    // find highest correlation coef.
    this->stat = calculateStat(this->genotype, pheno, &this->selected);
    this->perm.init(this->stat);

    // permutation
    double s;
    std::set<int> permSelected;
    while (this->perm.next()) {
      this->genotype.Transpose(genotype);
      permute(&pheno);

      s = calculateStat(this->genotype, pheno, &permSelected);
      this->perm.add(s);
    };
    fitOK = true;
    return (fitOK ? 0 : -1);
  };
  void reset() {
    this->perm.reset();
  }
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (isBinaryOutcome()){
      fprintf(fp, "NumIncMarker\t");
      this->perm.writeHeader(fp);
      fprintf(fp, "\n");
    } else
      fprintf(fp, "NA\n");
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (isBinaryOutcome()) {
      if (fitOK){
        fprintf(fp, "%zu\t", this->selected.size());
        this->perm.writeOutput(fp);
        fprintf(fp, "\n");
      } else {
        fprintf(fp, "NA\t");
        this->perm.writeOutput(fp);
        fprintf(fp, "\n");
      }
    } else {
      fprintf(fp, "NA\n");
    }
  };
  /**
   * For a given genotype and phenotype, calculate RareCover stats, which markers are selected
   * Here the genotype is: marker by people
   */
  double calculateStat(Matrix& genotype, Vector& phenotype, std::set<int>* selectedIndex) {
    std::set<int>& selected = *selectedIndex;
    selected.clear();

    Vector c; // collapsed genotype
    c.Dimension(phenotype.Length());
    c.Zero();

    double stat;
    while (selected.size() < genotype.rows) {
      int maxIdx = -1;
      double maxCorr = -1.0;
      double corr;
      for (int i = 0; i < genotype.rows; ++i){
        if (selected.count(i)) continue;
        corr = calculateCorrelation(genotype[i], c, phenotype);
        if (corr > maxCorr) {
          maxCorr = corr;
          maxIdx = i;
        }
      }
      if (maxIdx < 0) { // finish selection
        break;
      } else {
        if (maxCorr > stat) {
          // update selection
          stat = maxCorr;
          selected.insert(maxIdx);
          combine(&c, genotype[maxIdx]);
        } else { // no select any new marker
          break;
        };
      };
    }
    return stat;
  };
  /**
   * Calculate correlatio of (g + collapsed, pheno)
   */
  double calculateCorrelation(Vector& g, Vector& collapsed, Vector& pheno){
    double sum_g = 0.0;
    double sum_g2 = 0.0;
    double sum_p = 0.0;
    double sum_p2 = 0.0;
    double sum_gp = 0.0;
    int n = pheno.Length();
    for (int i = 0; i < n; ++i){
      double geno = (g[i] + collapsed[i] > 0 ) ? 1.0 : 0.0 ;
      if (geno > 0.0){
        sum_g += geno;
        sum_g2 += geno*geno;
        sum_gp += geno*pheno[i];
        sum_p += pheno[i];
        sum_p2 += pheno[i] * pheno[i];
      } else {
        sum_p += pheno[i];
        sum_p2 += pheno[i] * pheno[i];
      };
    };

    double cov_gp = sum_gp - sum_g * sum_p /n;
    double var_g = sum_g2 - sum_g * sum_g /n ;
    double var_p = sum_p2 - sum_p * sum_p /n;
    double v = var_g * var_p;
    if ( v < 1e-10) return 0.0;
    double corr = cov_gp / sqrt(v);
    return corr;
  };
  void combine(Vector* c, Vector& v){
    int n = v.Length();
    for (int i = 0; i < n; ++i){
      if ( (*c)[i] + v[i] > 0) {
        (*c)[i] = 1.0;
      }
    };
  };
private:
  Matrix genotype;
  std::set<int> selected;
  bool fitOK;
  int numVariant;
  double stat;
  Permutation perm;
}; //RareCoverTest

class CMATTest:public ModelFitter{
public:
CMATTest(int nPerm, double alpha): perm(nPerm, alpha) {
    this->modelName = "CMAT";
  }
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& covariate) {
    if (!isBinaryOutcome()) {
      fitOK = false;
      return -1;
    }
    if (covariate.cols = 0){
      fitOK = false;
      return -1;
    }
    copyPhenotype(phenotype, &this->pheno);

    // we use equal weight
    this->stat = this->calculateStat(genotype, pheno, &N_A, &N_U, &m_A, &m_U, &M_A, &M_U);
    this->perm.init(this->stat);

    // permutation part
    double d1,d2,d3,d4,d5,d6; // just used in permutation
    while(this->perm.next()) {
      permute(&pheno);
      // record new stats
      double pStat = this->calculateStat(genotype, pheno, &d1, &d2, &d3, &d4, &d5, &d6);
      this->perm.add(pStat);
    } // end permutation
    fitOK = true;
    return (fitOK ? 0 : -1);
  };
  void reset() {
    this->perm.reset();
  }
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (isBinaryOutcome()) {
      this->perm.writeHeader(fp);
      fprintf(fp, "\n");      
    } else
      fprintf(fp, "NA\n");
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (isBinaryOutcome()) {
      if (fitOK){
        this->perm.writeOutput(fp);
        fputs("\n", fp);
      } else {
        this->perm.writeOutput(fp);
        fputs("\n", fp);
      }
    } else {
      fprintf(fp, "\n");
    }
  };
  double calculateStat(Matrix& genotype, Vector& phenotype,
                       double* p_N_A, double* p_N_U,
                       double* p_m_A, double* p_m_U,
                       double* p_M_A, double* p_M_U){
    double& N_A = *p_N_A;
    double& N_U = *p_N_U;
    double& m_A = *p_m_A;
    double& m_U = *p_m_U;
    double& M_A = *p_M_A;
    double& M_U = *p_M_U;

    for (int i = 0; i < phenotype.Length(); ++i){
      if (phenotype[i] == 1) {
        ++ N_A;
      } else {
        ++ N_U;
      }
    }
    for (int i = 0; i < genotype.cols; ++i) {
      // for each marker, get its allele frequency
      double af = getMarkerFrequency(genotype, i);
      bool flip = false;
      if (af > 0.5) {
        flip = true;
      };
      for (int j = 0; j < genotype.rows; ++j) {
        if (phenotype[j] == 1) {
          if (!flip) {
            m_A += genotype[j][i];
            M_A += (2.0 - genotype[j][i]);
          } else {
            m_A += (2.0 - genotype[j][i]);
            M_A += genotype[j][i];
          }
        } else {
          if (!flip) {
            m_U += genotype[j][i];
            M_U += (2.0 - genotype[j][i]);
          } else {
            m_U += (2.0 - genotype[j][i]);
            M_U += genotype[j][i];
          }
        }
      };
    };
    if (N_A + N_U == 0.0) return 0.0;
    int numMarker = genotype.cols;
    return (N_A + N_U)/(2*N_A*N_U* numMarker) * (m_A*M_U - m_U*M_A)*(m_A*M_U - m_U*M_A)/(m_A+m_U)/(M_A+M_U);
  }
private:
  double N_A;
  double N_U;
  double m_A;
  double m_U;
  double M_A;
  double M_U;

  Vector pheno;
  bool fitOK;
  double stat;
  Permutation perm;
};//CMATTest

#if 0
class UStatTest{
};// UStatTest
#endif

/**
 * Implementation of Alkes Price's VT
 */
class VariableThresholdPrice: public ModelFitter{
public:
VariableThresholdPrice(int nPerm, double alpha): perm(nPerm, alpha) {
    this->modelName = "VariableThresholdPrice";
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& covariate) {
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    if (covariate.cols == 0) {
      fitOK = false;
      return -1;
    }

    rearrangeGenotypeByFrequency(genotype, &sortedGenotype, &this->freq);
    convertToReferenceAlleleCount(&sortedGenotype);
    collapseGenotype(&sortedGenotype);
    transpose(&sortedGenotype); // now each row is a collapsed genoype at certain frequency cutoff
    copyPhenotype(phenotype, &this->phenotype);

    this->zmax = -1.0;
    double z;
    if (isBinaryOutcome()) {
      for (int i = 0; i < sortedGenotype.rows; ++i) {
        z = calculateZForBinaryTrait(this->phenotype, this->sortedGenotype[i]);
        if ( z > this->zmax) {
          zmax = z;
          this->optimalFreq = freq[i];
        }
      }

      this->perm.init(zmax);

      // begin permutation
      while (this->perm.next()) {
        permute(&this->phenotype);
        double zp = -1.0;
        for (int j = 0; j < sortedGenotype.rows; ++j){
          z = calculateZForBinaryTrait(this->phenotype, this->sortedGenotype[j]);
          if (z > zp) {
            zp =z; 
          };
          if (zp > this->zmax)  // early stop
            break;
        }
        this->perm.add(zp);
      };

    } else {
      centerVector(&this->phenotype);
      for (int i = 0; i < sortedGenotype.rows; ++i) {
        z = calculateZForContinuousTrait(this->phenotype, this->sortedGenotype[i]);
        if ( z > this->zmax){
          this->zmax = z;
          this->optimalFreq = freq[i];
        }
      }
      this->perm.init(this->zmax);

      // begin permutation
      while(this->perm.next()) {
        double zp = -1.0;
        permute(&this->phenotype);
        for (int j = 0; j < sortedGenotype.rows; ++j){
          z = calculateZForContinuousTrait(this->phenotype, this->sortedGenotype[j]);
          if (z > zp) {
            zp = z;
          };
          if (zp > this->zmax) {
            break;
          }
        }
        this->perm.add(zp);
      };
    };

    fitOK = true;
    return 0;
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    this->perm.writeHeader(fp);
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    this->perm.writeOutput(fp);
    fprintf(fp, "\n");
  };
  void reset() {
    fitOK = true;
    this->perm.reset();
  };
private:
  /**
   * Convert genotype back to reference allele count
   * e.g. genotype 2 is reference allele count 0
   */
  void convertToReferenceAlleleCount(Matrix* g){
    Matrix& m = *g;
    for (int i = 0; i < m.rows; ++i) {
      for (int j = 0; j < m.cols; ++j) {
        m[i][j] = 2 - m[i][j];
      }
    }
  };
  /**
   * @param g is people by marker matrix
   * we will collpase left to right, column by column
   * it is mimic the behavior of setting different frequency cutoff
   */
  void collapseGenotype(Matrix* g){
    Matrix& m = *g;
    for (int i = 0; i < m.rows; ++i) {
      for (int j = 1; j < m.cols; ++j) {
        m[i][j] += m[i][j-1];
      }
    }
  }
  void transpose(Matrix* g) {
    Matrix tmp = *g;
    g->Transpose(tmp);
  };
  double calculateZForBinaryTrait(Vector& y, Vector& x){
    double ret = 0;
    int n = y.Length();
    for (int i = 0; i < n; ++i) {
      if (y[i]) ret += x[i];
    }
    return ret;
  };
  double calculateZForContinuousTrait(Vector& y, Vector& x){
    double ret = 0;
    int n = y.Length();
    for (int i = 0; i < n; ++i) {
      ret += x[i] * y[i];
    }
    return ret;
  };
private:
  Matrix sortedGenotype;
  std::vector<double> freq;
  bool fitOK;
  Vector phenotype;
  double zmax;
  double optimalFreq; // the frequency cutoff in unpermutated data which give smallest pvalue
  Permutation perm;
}; // VariableThresholdPrice


/**
 * This is variable threshold applying CMC test
 * NOTE: p-value is analytical (No permutation)
 *       Be cautious about the type-1 error
 */
class VariableThresholdCMC: public ModelFitter{
public:
VariableThresholdCMC():model(NULL),modelLen(0),modelCapacity(0){
    this->modelName = "VariableThresholdCMC";
    this->resize(32);
  };
  void resize(int n) {
    if (n < modelCapacity) {
      this->modelLen = n;
      return;
    }

    if (!model)
      delete[] model;

    model =  new CMCTest[n];
    this->modelLen = n;
    this->modelCapacity = n;
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& covariate) {
    if (genotype.cols > modelLen) {
      resize(genotype.cols);
      reset();
    }
    rearrangeGenotypeByFrequency(genotype, &sortedGenotype, &this->freq);
    for (int i = genotype.cols - 1; i >=0; --i){
      sortedGenotype.Dimension( genotype.rows, i + 1);
      if ( model[i].fit(phenotype, sortedGenotype, covariate) ) {
        fitOK = false;
      }
    }

    return 0;
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    model[0].writeHeader(fp, "FreqThreshold\t");
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    char buf[1000];
    // fputs(prependString, fp);
    for (int i = 0; i < freq.size(); i ++ ){
      sprintf(buf, "%s\t%f\t", prependString, freq[i]);
      model[i].writeOutput(fp, buf);
    }
  };
  void reset() {
    fitOK = true;
    for (int i = 0; i < this->modelLen; i++)
      model[i].reset();
  };
private:
  Matrix sortedGenotype;
  std::vector<double> freq;
  bool fitOK;
  CMCTest* model;
  int modelLen;
  int modelCapacity;
}; // VariableThresholdCMCTest

#if 0
class VariableThresholdFreqTest: public ModelFitter{
public:
  // write result header
  void writeHeader(FILE* fp) {
    fprintf(fp, "VT.FreqCutoff\tVT.PermPvalue");
  };
  // fitting model
  int fit(VCFData& data) {
    // arrange frequency
    for (int i = 0; i < data.collapsedMarkerFreq.size(); i++ ){
      order[ data.collapsedMarkerFreq[i] ] = i;
    }

    // collapsing using CMC
    progressiveMadsonBrowningCollapse(&data, &g, -1); // clear matrix

    // fit null model
    if (data.covariate && data.covariate->cols > 0) {
      fitOK = lrst.FitNullModel( *data.covariate, *data.extractPhenotype(), 100);
      if (!fitOK)
        return -1;
    }

    for ( order_it = order.begin();
          order_it != order.end();
          order_it++ ){
      progressiveMadsonBrowningCollapse(&data, &g, order_it->second);
      outFreq.push_back(order_it->first);

      if (data.covariate && data.covariate->cols > 0) {
        fitOK = lrst.TestCovariate( *data.covariate, *data.extractPhenotype(), g);
        if (!fitOK) {
          pvalue.push_back(-1.0);
          continue;
        }
      } else {
        fitOK = lrst.TestCovariate( *data.collapsedGenotype, *data.extractPhenotype());
        if (!fitOK) {
          pvalue.push_back(-1.0);
          continue;
        }
      }
      pvalue.push_back(lrst.GetPvalue());
    }

    return 0;
  };
  // write model output
  void writeOutput(FILE* fp) {
    assert(outFreq.size() == pvalue.size());

    for (int i = 0; i < outFreq.size(); i++) {
      if (i) fputc(',', fp);
      fprintf(fp, "%g", outFreq[i]);
    }
    fputc('\t', fp);
    for (int i = 0; i < pvalue.size(); i++) {
      if (i) fputc(',', fp);
      if (pvalue[i] > 0)
        fprintf(fp, "%g", pvalue[i]);
      else
        fputs("NA", fp);
    }
  };
  void reset() {
    this->outFreq.clear();
    this->pvalue.clear();
    this->order.clear();
  };
private:
  Matrix g;
  LogisticRegressionScoreTest lrst;
  std::vector<double> outFreq;
  std::vector<double> pvalue;
  std::map<double, int> order;
  std::map<double, int>::const_iterator order_it;
  bool fitOK;
}; // VariableThresholdFreqTest
#endif

class SkatTest: public ModelFitter{
public:
  /* SkatTest(const std::vector<std::string>& param) { */
SkatTest(int nPerm, int alpha, double beta1, double beta2):perm(nPerm, alpha) {
    if (nPerm >0)
      this->usePermutation = true;
    this->beta1 = beta1;
    this->beta2 = beta2;
    this->modelName = "Skat";
  };
  void reset() {
    this->skat.Reset();
    this->perm.reset();
    stat = -9999;
  }
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& covariate) {
    // fill it wegith
    // NOTE: our frequency calculation is slightly different than SKAT, so we will need to adjust it back
    weight.Dimension(genotype.cols);
    for (int i = 0; i < weight.Length(); i++) {
      double freq = getMarkerFrequency(genotype, i);
      freq = (freq * (genotype.rows * 2 + 1) -1 ) / (genotype.rows *2 );
      if (freq > 1e-30) { // avoid dividing zero
        weight[i] =  gsl_ran_beta_pdf(freq, this->beta1, this->beta2);  /// default SKAT use beta(MAF, 1, 25)
        weight[i] *= weight[i];
        // fprintf(stderr, "weight(%d, %d, %f ) = %f\n", 1, 25, freq, weight[i]);
      } else {
        weight[i] = 0.0;
      }
    };

    // get ynulll
    Vector phenoVec;
    copyPhenotype(phenotype, &phenoVec);

    // ynull is mean of y (removing genotypes) in model Ynull ~ X (aka Ynull ~ X + 0.0 * G )
    Matrix cov;
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);

    if (isBinaryOutcome()) {
      fitOK = logistic.FitLogisticModel(cov, phenoVec, 100);
      if (!fitOK) {
        return -1;
      }
      ynull = logistic.GetPredicted();
      v = logistic.GetVariance();
    } else {
      fitOK = linear.FitLinearModel(cov, phenoVec);
      if (!fitOK) {
        return -1;
      }
      ynull = linear.GetPredicted();
      v.Dimension(genotype.rows);
      for (int i = 0; i < genotype.rows; ++i) {
        v[i] = linear.GetSigma2();
      }
    }

    // get Pvalue
    skat.CalculatePValue(phenoVec, ynull, cov, v, genotype, weight);
    this->pValue = skat.GetPvalue();

    // permuation part
    this->stat =  skat.GetQ();
    this->perm.init(this->stat);

    double s;
    while (this->perm.next()) {
      permute(&phenoVec);
      skat.CalculatePValue(phenoVec, ynull, cov, v, genotype, weight);
      s = skat.GetQ();
      this->perm.add(s);
    };
    fitOK = true;
    return 0;
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (!usePermutation)
      fprintf(fp, "NMARKER\tQ\tPvalue\n");
    else {
      fprintf(fp, "NMARKER\tQ\tPvalue\t");
      this->perm.writeHeader(fp);
      fprintf(fp, "\n");
    }
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (!fitOK){
      fprintf(fp, "%d\tNA\tNA\n", this->weight.Length());
    } else {
      // binary outcome and quantative trait are similar output
      fprintf(fp, "%d\t%g\t%g", this->weight.Length(), this->skat.GetQ(), this->pValue);
      if (usePermutation) {
        fprintf(fp, "\t");
        this->perm.writeOutput(fp);
      }
      fprintf(fp, "\n");
    }
  };
private:
  double beta1;
  double beta2;
  // Matrix X; // n by (p+1) matrix, people by covariate (note intercept is needed);
  Vector v;
  Vector weight;
  LogisticRegression logistic;
  LinearRegression linear;
  Vector ynull;
  Skat skat;
  bool fitOK;
  double pValue;
  int nMarker;

  bool usePermutation;
  double stat;
  Permutation perm;
}; // SkatTest

class KbacTest: public ModelFitter{
public:
KbacTest(int nPerm):nPerm(nPerm), xdatIn(NULL), ydatIn(NULL),mafIn(NULL), xcol(0), ylen(0), nn(0), qq(0) {
    this->modelName = "KBAC";
  };
  ~KbacTest() {
    if (this->xdatIn) delete[] this->xdatIn;
    if (this->ydatIn) delete[] this->ydatIn;
    if (this->mafIn) delete[] this->mafIn;
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fprintf(fp, "KBAC.Pvalue\n");
  };
  void reset() {
    // clear_kbac_test();
  }
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& covariate) {
    if (!isBinaryOutcome()) {
      fitOK = false;
      return -1;
    }
    if (covariate.cols  != 0){
      fitOK = false;
      return -1;
    }

    this->resize(genotype.rows, genotype.cols);
    this->nn = this->nPerm;
    this->qq = 1;
    this->aa = 0.05;
    this->mafUpper = 1.0; // no need to further prune alleles
    this->twosided = 1;
    // genotype is: people by marker
    for (int i = 0; i < genotype.rows; ++i) {
      for (int j = 0; j < genotype.cols; ++j) {
        // note: KBAC package main.R
        // we essentially store xmat in a array where the order is
        // people1 marker1, people1 marker2, people1 marker3.... then
        // people2 marker1, people2 marker2, people2 marker3....
        xdatIn[i * genotype.cols + j] = genotype[i][j];
        /* if (genotype[i][j] != 0.0) { */
        /*   fprintf(stderr, "i=%d, j=%d, genotype=%g\n", i,j,genotype[i][j]); */
        /*   fprintf(stderr, "j*genotype.rows +i = %d, xdatIn[..] = %g\n", j * genotype.rows + i, xdatIn[j * genotype.rows + i]); */
        /* } */
      }
    }
    for (int i = 0; i < genotype.rows; ++i ) {
      ydatIn[i] = phenotype[i][0];
    }
    for (int j = 0; j < genotype.cols; ++j ) {
      mafIn[j] = getMarkerFrequency(genotype, j);
    }

    /**
     * nn: number of permutation
     * qq: quiet
     * aa: alpha level
     * mafUpper: MAF upper threshold
     * xdatIn: genotype matrix
     * ydatIn: phenotype matrix
     * mafIn:  allele frequency for each marker of x
     * xcol: number of column of genotype
     * ylen: length of y
     */
    set_up_kbac_test(&nn, &qq, &aa, &mafUpper, xdatIn, ydatIn, mafIn,  &xcol, &ylen);
    /**
     * pvalue: results are here
     * twosided: two sided tests or one sided test
     */
    this->pValue = 9.0; // it is required to pass pvalue = 9.0 (see interface.R)
    do_kbac_test(&this->pValue, &this->twosided);
    // fprintf(stderr, "pValue = %g, twoside = %d\n", this->pValue, this->twosided);
    clear_kbac_test();
    this->fitOK = true;
    return 0;
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (!fitOK){
      fputs("NA\n", fp);
    } else {
      if (isBinaryOutcome() ) {
        fprintf(fp, "%f\n", this->pValue);
      } else {
        fprintf(fp, "%f\n", this->pValue);
      }
    }
  };
  void resize(int numPeople, int numMarker) {
    bool resized = false;
    if (numPeople != this->ylen) {
      delete[] this->ydatIn;
      this->ydatIn = new double[ numPeople];
      this->ylen = numPeople;
      resized = true;
    }
    if (numMarker != this->xcol) {
      delete[] this->mafIn;
      this->mafIn = new double[ numMarker];
      this->xcol = numMarker;
      resized = true;
    }

    if (resized) {
      delete[] this->xdatIn;
      this->xdatIn = new double[ numPeople * numMarker ];
    }
  };
private:
  int nPerm;

  int nn;
  int qq;
  double aa;
  double mafUpper;
  double* xdatIn;
  double* ydatIn;
  double* mafIn;
  int xcol;
  int ylen;

  int twosided;
  bool fitOK;
  double pValue;
}; // KbacTest

class DumpModel: public ModelFitter{
public:
  DumpModel(const char* prefix) {
    this->prefix = prefix;
    this->modelName = "DumpData";
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) { // e.g. column headers.
    fputs(prependString, fp);
    fprintf(fp, "FileName");

    this->header = prependString;
  }
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype, Matrix& covariate) {
    this->phenotype = phenotype;
    this->genotype = genotype;
    this->covariate = covariate;
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString){
    std::string fn = this->prefix + "\t" + prependString + "data";
    for (int i = 0; i < fn.size(); ++i) {
      if (fn[i] == '\t') fn[i] = '.';
    }

    fputs(prependString, fp);
    fprintf(fp, "%s\n", fn.c_str());

    // write header
    FILE* fDump= fopen(fn.c_str(), "wt");
    fprintf(fDump, "%s\t", this->header.c_str());
    if (phenotype.cols == 1)  {
      fprintf(fDump, "Y");
    } else {
      for (int i = 0; i < phenotype.cols; i++) {
        if (i) fprintf(fDump, "\t");
        fprintf(fDump, "Y%d", i);
      }
    }
    for (int i = 0; i < genotype.cols; i++) {
      fprintf(fDump, "\tX%d", i);
    };
    for (int i = 0; i < covariate.cols; i++) {
      fprintf(fDump, "\tC%d", i);
    };
    fprintf(fDump, "\n");

    // write content
    for (int i = 0; i < phenotype.rows; ++i) {
      fputs(prependString, fDump);
      for (int j = 0; j < phenotype.cols; ++j) {
        if (j) fprintf(fDump, "\t");
        fprintf(fDump, "%f", phenotype[i][j]);
      }
      for (int j = 0; j < genotype.cols; ++j) {
        fprintf(fDump, "\t%f", genotype[i][j]);
      }
      for (int j = 0; j < covariate.cols; ++j) {
        fprintf(fDump, "\t%f", covariate[i][j]);
      }
      fprintf(fDump, "\n");
    };
    fclose(fDump);

  };

  void reset() {
  }; // for particular class to call when fitting repeatedly
private:
  Matrix phenotype;
  Matrix genotype;
  Matrix covariate;
  std::string prefix;
  std::string header;
}; // end DumpModel



//////////////////////////////////////////////////////////////////////
// Implementation of various collpasing methods
/**
 * @return Madson-Browning definition of alleleFrequency
 */
double getMarkerFrequency(Matrix& in, int col){
  int& numPeople = in.rows;
  int ac = 0;
  int an = 0;
  for (int p = 0; p < numPeople; p++) {
    if (in[p][col] >= 0) {
      ac += in[p][col];
      an += 2;
    }
  }
  if ( an == 0 ) return 0.0;
  //double freq = 1.0 * (ac + 1) / (an + 1);
  double freq = 1.0 * ac / an;
  return freq;
};

double getMarkerFrequencyFromControl(Matrix& in, Vector& pheno, int col){
  int& numPeople = in.rows;
  int ac = 0;
  int an = 0;
  for (int p = 0; p < numPeople; p++) {
    if (pheno[p] == 1) continue;
    if (in[p][col] >= 0) {
      ac += in[p][col];
      an += 2;
    }
  }
  double freq = 1.0 * (ac + 1) / (an + 1);
  return freq;
};

void cmcCollapse(Matrix& in, Matrix* out){
  assert(out);
  int numPeople = in.rows;
  int numMarker = in.cols;

  out->Dimension(numPeople, 1);
  out->Zero();
  for (int p = 0; p < numPeople; p++){
    for (int m = 0; m < numMarker; m++) {
      int g = (int)(in[p][m]);
      if (g > 0) {
        (*out)[p][0] = 1.0;
        break;
      }
    };
  };
};

void zegginiCollapse(Matrix& in, Matrix* out){
  assert(out);
  int numPeople = in.rows;
  int numMarker = in.cols;

  out->Dimension(numPeople, 1);
  out->Zero();
  for (int p = 0; p < numPeople; p++){
    for (int m = 0; m < numMarker; m++) {
      int g = (int)(in[p][m]);
      if (g > 0) { // genotype is non-reference
        (*out)[p][0] += 1.0;
      }
    };
  };
};

/**
 * @param genotype : people by marker matrix
 * @param phenotype: binary trait (0 or 1)
 * @param out: collapsed genotype
 */
void madsonBrowningCollapse(Matrix& genotype, Vector& phenotype, Matrix* out){
  assert(out);
  int& numPeople = genotype.rows;
  int numMarker = genotype.cols;

  out->Dimension(numPeople, 1);
  out->Zero();

  for (int m = 0; m < numMarker; m++) {
    // calculate weight
    double freq = getMarkerFrequencyFromControl(genotype, phenotype, m);
    if (freq <= 0.0 || freq >= 1.0) continue; // avoid freq == 1.0
    double weight = 1.0 / sqrt(freq * (1.0-freq));
    // fprintf(stderr, "freq = %f\n", freq);

    for (int p = 0; p < numPeople; p++) {
      (*out)[p][0] += genotype[p][m] * weight;
    }
  };
};

void fpCollapse(Matrix& in, Matrix* out){
  assert(out);
  int& numPeople = in.rows;
  int numMarker = in.cols;

  out->Dimension(numPeople, 1);
  out->Zero();

  for (int m = 0; m < numMarker; m++) {
    // calculate weight
    double freq = getMarkerFrequency(in, m);
    if (freq <= 0.0 || freq >= 1.0) continue; // avoid freq == 1.0
    double weight = 1.0 / sqrt(freq * (1.0-freq));
    // fprintf(stderr, "freq = %f\n", freq);

    for (int p = 0; p < numPeople; p++) {
      (*out)[p][0] += in[p][m] * weight;
    }
  };
};

void madsonBrowningCollapse(Matrix* d, Matrix* out){
  assert(out);
  Matrix& in = (*d);
  int& numPeople = in.rows;
  int numMarker = in.cols;

  out->Dimension(numPeople, 1);
  out->Zero();


  for (int m = 0; m < numMarker; m++) {
    // calculate weight
    double freq = getMarkerFrequency(in, m);
    if (freq <= 0.0 || freq >= 1.0) continue; // avoid freq == 1.0
    double weight = 1.0 / sqrt(freq * (1.0-freq));
    // fprintf(stderr, "freq = %f\n", freq);

    for (int p = 0; p < numPeople; p++) {
      (*out)[p][0] += in[p][m] * weight;
    }
  };
};

template<class T>
class OrderFunction {
public:
OrderFunction(T& t): v(t) {};
  bool operator() (int i, int j)  {
    return v[i] < v[j];
  };
  const T& v;
};

/**
 * @param freq: 0.3, 0.2, 0.1, 0.4
 * will return
 * @param ord:  2, 1, 0, 3
 *
 * algorithm:
 * make a pair like this: (0.3, 0), (0.2, 1), (0.1, 2), (0.4, 3)
 * sort by first element:
 * (0.1, 2), (0.2, 1), (0.3, 0), (0.4, 3)
 * extract second element:
 * 2, 1, 0, 3
 */
void order(std::vector<double>& freq, std::vector<int>* ord) {
  ord->resize(freq.size());
  for (int i = 0; i < freq.size(); ++i)
    (*ord)[i] = i;

  OrderFunction< std::vector<double> > func(freq);
  std::sort(ord->begin(), ord->end(), func);
};

/**
 * Collapsing @param d to @param out, the order of columns in @param out is the same as @param freq
 * which is the frequency upper bounds, and @param freq will be increase frequency.
 * e.g.
 * @param in P by 3 matrix, then @param out will be 3 columns too
 * if @param freq = (0.1, 0.2, 0.3) meaning
 * @param out column 0 using 0.1 frequency upper bound, and the smallest frequency of marker is 0.1 (@param freq[0])
 * @param out column 1 using 0.2 frequency upper bound, and the second largest frequency of marker is 0.2 (@param freq[1])
 * @param out column 2 using 0.3 frequency upper bound, and the largest frequency of marker is 0.3 (@param freq[2])
 */
void rearrangeGenotypeByFrequency(Matrix& in, Matrix* out, std::vector<double>* freq) {
  assert(out && freq);
  out->Dimension(in.rows, in.cols);

  const int& numPeople = in.rows;
  const int& numMarker = in.cols;
  freq->resize(numMarker);

  /* out->Dimension(numPeople, numMarker); */
  /* if (col < 0) { */
  /*   out->Zero(); */
  /*   return; */
  /* } */

  for (int m = 0; m < numMarker; ++m){
    (*freq)[m] = getMarkerFrequency(in, m);
  }

  std::vector<int> ord;
  order(*freq, &ord);
  std::sort(freq->begin(), freq->end());

  for (int m = 0; m < numMarker; ++m){
    const int& col = ord[m];
    for (int p = 0; p < numPeople; ++p) {
      (*out)[p][m] = in[p][col];
    }
  }
};

void progressiveCMCCollapse(Matrix* d, Matrix* out, std::vector<double>* freq) {
  assert(out && freq);
  Matrix& in = (*d);
  out->Dimension(in.rows, in.cols);

  const int& numPeople = in.rows;
  const int& numMarker = in.cols;
  freq->resize(numMarker);

  /* out->Dimension(numPeople, numMarker); */
  /* if (col < 0) { */
  /*   out->Zero(); */
  /*   return; */
  /* } */

  for (int m = 0; m < numMarker; ++m){
    (*freq)[m] = getMarkerFrequency(in, m);
  }

  std::vector<int> ord;
  order(*freq, &ord);
  std::sort(freq->begin(), freq->end());

  for (int m = 0; m < numMarker; ++m){
    const int& col = ord[m];
    if (m == 0) {
      for (int p = 0; p < numPeople; ++p) {
        if (in[p][col] > 0) {
          (*out)[p][m] = 1;
        }
      }
    } else { //
      for (int p = 0; p < numPeople; ++p) {
        if ((*out)[p][m-1] > 0) {
          (*out)[p][m] = 1;
          continue;
        };
        if (in[p][col] > 0) {
          (*out)[p][m] = 1;
          continue;
        }
      }
    }
  }

};

void progressiveMadsonBrowningCollapse(Matrix* d, Matrix* out, std::vector<double>* freq) {
  assert(out);
  Matrix& in = (*d);
  out->Dimension(in.rows, in.cols);

  int numPeople = in.rows;
  int numMarker = in.cols;
  freq->resize(numMarker);

  /* out->Dimension(numPeople, 1); */
  /* if (col < 0) { */
  /*   out->Zero(); */
  /*   return; */
  /* } */

  for (int m = 0; m < numMarker; ++m){
    (*freq)[m] = getMarkerFrequency(in, m);
  }

  std::vector<int> ord;
  order(*freq, &ord);
  std::sort(freq->begin(), freq->end());

  double weight;
  for (int m = 0; m < numMarker; ++m){
    const int& col = ord[m];
    weight = 1.0 / sqrt (  (*freq)[m] * (1.0 - (*freq)[m])) ;

    if (m == 0) {
      for (int p = 0; p < numPeople; ++p) {
        (*out)[p][m] = weight * in[p][col];
      }
    } else { //
      for (int p = 0; p < numPeople; ++p) {
        (*out)[p][m] = (*out)[p][m-1] + weight * in[p][col];
      }
    }
  }
};

void permute(Vector* vec){
  Vector& v = *vec;
  int n = v.Length();
  double tmp;
  for (int i = n - 1; i >= 1; --i) {
    // pick j from 0 <= j <= i
    int j = rand() % (i+1);
    if (i != j) {
      tmp = v[i];
      v[i] = v[j];
      v[j] = tmp;
    }
  }
};

void centerVector(Vector* v){
  double avg = v->Average();
  int n = v->Length();
  for (int i = 0; i < n; ++i) {
    (*v)[i] -= avg;
  };
};


#endif /* _MODELFITTER_H_ */
