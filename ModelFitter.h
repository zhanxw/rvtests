#ifndef _MODELFITTER_H_
#define _MODELFITTER_H_

#include "libsrc/MathMatrix.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "regression/LogisticRegression.h"
#include "regression/LogisticRegressionScoreTest.h"
#include "regression/LinearRegression.h"
#include "regression/LinearRegressionScoreTest.h"
#include "regression/Skat.h"
#include "regression/Table2by2.h"
#include "regression/kbac_interface.h"
#include "regression/FastLMM.h"

#include "regression/MatrixOperation.h"
#include "DataConsolidator.h"
#include "LinearAlgebra.h"
#include "snp_hwe.c"
#include "Result.h"
#include "Summary.h"
#include "Permutation.h"

extern SummaryHeader* g_SummaryHeader;

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

// take X, Y, Cov and fit model
// note, ModelFitter will use VCFData as READ-ONLY data structure,
// and collapsing results are stored internally.
class ModelFitter{
public:
  virtual int fit(DataConsolidator* dc) = 0;

  // write result header
  virtual void writeHeader(FileWriter* fp, const Result& siteInfo) = 0;
  // write model output
  virtual void writeOutput(FileWriter* fp, const Result& siteInfo) = 0;

  ModelFitter(){
    this->modelName = "Unassigned_Model_Name";
    this->binaryOutcome = false; // default: using continuous outcome
  };
  virtual ~ModelFitter() {};
  const std::string& getModelName() const { return this->modelName; };
  // for particular class to call when fitting repeatedly
  // e.g. clear permutation counter
  // e.g. clear internal cache
  virtual void reset() { this->result.clearValue();};
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
  const Result& getResult() const {
    return this->result;
  };
protected:
  std::string modelName;
  bool binaryOutcome;
  Result result;
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
    result.addHeader("Test");
    result.addHeader("Beta");
    result.addHeader("SE");
    result.addHeader("Pvalue");

  };
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& cov= dc->getCovariate();

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
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    // fprintf(fp, "Test\tBeta\tSE\tPvalue\n");
    result.writeHeaderLine(fp);
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    // skip interecept (column 0)
    for (int i = 1; i < this->X.cols; ++i) {
      siteInfo.writeValueTab(fp);
      if (!fitOK) {
        // fprintf(fp, "%s\tNA\tNA\tNA\n", this->X.GetColumnLabel(i));
        result.updateValue("Test", this->X.GetColumnLabel(i));
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

        // fprintf(fp, "%s\t%g\t%g\t%g\n", this->X.GetColumnLabel(i), beta, se, pval);
        result.updateValue("Test", this->X.GetColumnLabel(i));
        result.updateValue("Beta", beta);
        result.updateValue("SE", se);
        result.updateValue("Pvalue", pval);
      }
      result.writeValueLine(fp);
    }
  };
private:
  Matrix X; // 1 + cov + geno
  Vector Y; // phenotype
  LinearRegression linear;
  LogisticRegression logistic;
  bool fitOK;
}; // SingleVariantWaldTest

/* /\** */
/*  * @return 0 if success */
/*  *\/ */
/* int checkHWEandCallRate(Matrix& geno, int* homRef, int* het, int* homAlt, double* hweP, double* callRate) { */
/*   *homRef = *het = *homAlt = 0; */
/*   *hweP = 0.0; */
/*   *callRate = 0.0; */

/*   int nonMissingGeno = 0; */
/*   int totalGeno = geno.rows; */

/*   for (int i = 0; i < totalGeno; ++i) { */
/*     const int& g = geno[i][0]; */
/*     if (g < 0) continue; */
/*     if (g == 0) { */
/*       ++ (*homRef); */
/*       ++ (nonMissingGeno); */
/*     } else if (g == 1) { */
/*       ++ (*het); */
/*       ++ (nonMissingGeno); */
/*     } else if (g == 2) { */
/*       ++ (*homAlt); */
/*       ++ (nonMissingGeno); */
/*     } */
/*   } */
/*   if (totalGeno) */
/*     (*callRate) = 1.0 * nonMissingGeno / totalGeno; */
/*   if ( (*homRef) >= 0 && (*het) >= 0 && (*homAlt) >= 0) */
/*     (*hweP) = SNPHWE( (*het), (*homRef), (*homAlt)); */
/*   return 0; */
/* } */

class SingleVariantScoreTest: public ModelFitter{
public:
  SingleVariantScoreTest(){
    this->modelName = "SingleScore";
  };
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate= dc->getCovariate();

    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }
    nSample = genotype.rows;
    af = getMarkerFrequency(genotype, 0);

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
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    if (g_SummaryHeader) {
      g_SummaryHeader->outputHeader(fp);
    }

    siteInfo.writeHeaderTab(fp);
    // fprintf(fp, "AF\tStat\tDirection\tPvalue\n");
    result.addHeader("AF");
    result.addHeader("STAT");
    result.addHeader("DIRECTION");
    result.addHeader("PVALUE");
    result.writeHeaderLine(fp);
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    result.clearValue();
    if (fitOK) {
      if (!isBinaryOutcome()) {
        result.updateValue("AF", af);
        result.updateValue("STAT", linear.GetStat());
        result.updateValue("DIRECTION", linear.GetU()[0][0] > 0 ? "+": "-");
        result.updateValue("PVALUE", linear.GetPvalue());
      } else {
        result.updateValue("AF", af);
        result.updateValue("STAT", logistic.GetStat());
        result.updateValue("DIRECTION", logistic.GetU()[0][0] > 0 ? "+" : "-");
        result.updateValue("PVALUE", logistic.GetPvalue());
      }
    }
    result.writeValueLine(fp);
  };
private:
  double af;
  int nSample;
  Vector pheno;
  LinearRegressionScoreTest linear;
  LogisticRegressionScoreTest logistic;
  bool fitOK;
  Matrix cov;
}; // SingleVariantScoreTest

class SingleVariantFisherExactTest: public ModelFitter{
public:
  SingleVariantFisherExactTest() {
    this->modelName = "FisherExact";
  };
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.addHeader("Fisher.N00");
    result.addHeader("Fisher.N01");
    result.addHeader("Fisher.N10");
    result.addHeader("Fisher.N11");
    result.addHeader("CtrlAF");
    result.addHeader("CaseAF");
    result.addHeader("Fisher.PvalueTwoSide");
    result.addHeader("Fisher.PvalueLess");
    result.addHeader("Fisher.PvalueGreater");
    result.writeHeaderLine(fp);
  };
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& cov= dc->getCovariate();

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
    return 0;
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (fitOK) {
      /* fprintf(fp, "%d\t", model.Get00()); */
      /* fprintf(fp, "%d\t", model.Get01()); */
      /* fprintf(fp, "%d\t", model.Get10()); */
      /* fprintf(fp, "%d\t", model.Get11()); */
      result.updateValue("Fisher.N00", model.Get00());
      result.updateValue("Fisher.N01", model.Get01());
      result.updateValue("Fisher.N10", model.Get10());
      result.updateValue("Fisher.N11", model.Get11());

      if (ctrlAN == 0) {
        // fprintf(fp, "0\t");
        result.updateValue("CtrlAF", 0);
      } else{
        // fprintf(fp, "%g\t", 1.0 * ctrlAC / ctrlAN);
        result.updateValue("CtrlAF", 1.0 * ctrlAC / ctrlAN);
      }
      if (caseAN == 0) {
        // fprintf(fp, "0\t");
        result.updateValue("CaseAF", 0);
      } else{
        // fprintf(fp, "%g\t", 1.0 * caseAC / caseAN);
        result.updateValue("CaseAF", 1.0 * caseAC / caseAN);
      }
      /* fprintf(fp, "%lf\t"  , model.getPExactTwoSided()); */
      /* fprintf(fp, "%lf\t"  , model.getPExactOneSidedLess()); */
      /* fprintf(fp, "%lf\n"  , model.getPExactOneSidedGreater()); OB*/
      result.updateValue("Fisher.PvalueTwoSide", model.getPExactTwoSided());
      result.updateValue("Fisher.PvalueLess", model.getPExactOneSidedLess());
      result.updateValue("Fisher.PvalueGreater", model.getPExactOneSidedGreater());
    } /* else { */
    /*   fprintf(fp, "0\t0\t0\t0\tNA\tNA\tNA\n"); */
    /* } */
    result.writeValueLine(fp);
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

class SingleVariantFamilyScore: public ModelFitter{
public:
SingleVariantFamilyScore():model(FastLMM::SCORE, FastLMM::MLE) {
    this->modelName = "FamScore";
    result.addHeader("AF");
    result.addHeader("U.Stat");
    result.addHeader("V.Stat");
    result.addHeader("FamScore.Pvalue");
    needToFitNullModel = true;
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    if (isBinaryOutcome()) {
      fitOK = false;
      return -1;
    }
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();

    if (needToFitNullModel || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK = (0 == model.FitNullModel(cov, phenotype, *dc->getKinshipU(), *dc->getKinshipS()) ? true: false);
      if (!fitOK) return -1;
      needToFitNullModel = false;
    }
    
    fitOK = (0 == model.TestCovariate(cov, phenotype, genotype, *dc->getKinshipU(), *dc->getKinshipS()) ? true: false);
    af = model.GetAF(*dc->getKinshipU(), *dc->getKinshipS());
    u = model.GetUStat();
    v = model.GetVStat();
    pvalue = model.GetPValue();
    return (fitOK ? 0 : -1);
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.writeHeaderLine(fp);
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (fitOK) {
      // result.updateValue("NonRefSite", this->totalNonRefSite());
      if (isBinaryOutcome()) {
        //result.updateValue("CMC.Pvalue", logistic.GetPvalue());
      } else {
        result.updateValue("AF", af);
        result.updateValue("U.Stat", u);
        result.updateValue("V.Stat", v);
        result.updateValue("FamScore.Pvalue", pvalue);
      }
    }
    result.writeValueLine(fp);
  };
  
private:
  Matrix cov;
  FastLMM model;
  bool needToFitNullModel;
  bool fitOK;
  double af;
  double u;
  double v;
  double pvalue;
};

class CMCTest: public ModelFitter{
public:
  CMCTest() {
    this->modelName = "CMC";
    result.addHeader("NonRefSite");
    result.addHeader("CMC.Pvalue");
  };
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();
    return this->fit(phenotype, genotype, covariate);
  }
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
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.writeHeaderLine(fp);
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (fitOK) {
      result.updateValue("NonRefSite", this->totalNonRefSite());
      if (isBinaryOutcome()) {
        result.updateValue("CMC.Pvalue", logistic.GetPvalue());
      } else {
        result.updateValue("CMC.Pvalue", linear.GetPvalue());
      }
    }
    result.writeValueLine(fp);
  };
private:
  /**
   * If the genotype is not exactly 0.0, we will count is as non-reference site
   */
  int totalNonRefSite(){
    int s = 0;
    for (int i = 0; i < collapsedGenotype.rows; ++i) {
      s += collapsedGenotype[i][0] == 0.0 ? 0 : 1;
    }
    return (s);
  }

  Matrix collapsedGenotype;
  LogisticRegressionScoreTest logistic;
  LinearRegressionScoreTest linear;
  bool fitOK;
  int numVariant;
}; // CMCTest

class CMCWaldTest: public ModelFitter{
public:
  CMCWaldTest() {
    this->modelName = "CMCWald";
    result.addHeader("NonRefSite");
    result.addHeader("Beta");
    result.addHeader("SE");
    result.addHeader("CMCWald.Pvalue");
  };
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();

    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }

    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i< phenotype.rows; i++){
      pheno[i] = phenotype[i][0];
    }

    cmcCollapse(genotype, &collapsedGenotype);

    if (covariate.cols) {
      copyGenotypeWithCovariateAndIntercept(collapsedGenotype, covariate, &this->X);
    } else {
      copyGenotypeWithIntercept(collapsedGenotype, &this->X);
    }

    if (!isBinaryOutcome()) {
      fitOK = linear.FitLinearModel(X, phenotype);
    } else {
      fitOK = logistic.FitLogisticModel(X, phenotype, 100);
    }
    return (fitOK ? 0 : -1);
  };
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.writeHeaderLine(fp);
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    for (int i = 1; i < this->X.cols; ++i) {
      siteInfo.writeValueTab(fp);
      if (fitOK) {
        result.updateValue("NonRefSite", this->totalNonRefSite());
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
        result.updateValue("Beta", beta);
        result.updateValue("SE", se);
        result.updateValue("CMCWald.Pvalue", pval);
      }
      result.writeValueLine(fp);
    }
  };
private:
  /**
   * If the genotype is not exactly 0.0, we will count is as non-reference site
   */
  int totalNonRefSite(){
    int s = 0;
    for (int i = 0; i < collapsedGenotype.rows; ++i) {
      s += collapsedGenotype[i][0] == 0.0 ? 0 : 1;
    }
    return (s);
  }

  Matrix collapsedGenotype;
  Matrix X;
  LogisticRegression logistic;
  LinearRegression linear;
  bool fitOK;
  int numVariant;
}; // CMCWaldTest

class CMCFisherExactTest: public ModelFitter{
public:
  CMCFisherExactTest() {
    this->modelName = "CMCFisherExact";
    result.addHeader("exactCMC.N00");
    result.addHeader("exactCMC.N01");
    result.addHeader("exactCMC.N10");
    result.addHeader("exactCMC.N11");
    result.addHeader("exactCMC.PvalueTwoSide");
    result.addHeader("exactCMC.PvalueLess");
    result.addHeader("exactCMC.PvalueGreater");
  };
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& cov= dc->getCovariate();

    if (!isBinaryOutcome()){
      fitOK = false;
      return -1;
    }
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
    return 0;
  };
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.writeHeaderLine(fp);
    /* fputs("exactCMC.N00\t", fp); */
    /* fputs("exactCMC.N01\t", fp); */
    /* fputs("exactCMC.N10\t", fp); */
    /* fputs("exactCMC.N11\t", fp); */
    /* fputs("exactCMC.PvalueTwoSide\t", fp); */
    /* fputs("exactCMC.PvalueLess\t", fp); */
    /* fputs("exactCMC.PvalueGreater\n", fp); */
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    result.clearValue();
    if (fitOK) {
      /* fprintf(fp, "%d\t", model.Get00()); */
      /* fprintf(fp, "%d\t", model.Get01()); */
      /* fprintf(fp, "%d\t", model.Get10()); */
      /* fprintf(fp, "%d\t", model.Get11()); */

      /* fprintf(fp, "%lf\t", model.getPExactTwoSided()); */
      /* fprintf(fp, "%lf\t", model.getPExactOneSidedLess()); */
      /* fprintf(fp, "%lf\n"  , model.getPExactOneSidedGreater()); */

      result.updateValue("exactCMC.N00", model.Get00());
      result.updateValue("exactCMC.N01", model.Get01());
      result.updateValue("exactCMC.N10", model.Get10());
      result.updateValue("exactCMC.N11", model.Get11());
      result.updateValue("exactCMC.PvalueTwoSide", model.getPExactTwoSided());
      result.updateValue("exactCMC.PvalueLess", model.getPExactOneSidedLess());
      result.updateValue("exactCMC.PvalueGreater", model.getPExactOneSidedGreater());

    }
    /*  else { */
    /*   fprintf(fp, "0\t0\t0\t0\tNA\tNA\tNA\n"); */
    /* } */
    result.writeValueLine(fp);
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
    this->modelName = "Zeggini";
    result.addHeader("Zeggini.Pvalue");
  };
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();

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
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    //fprintf(fp, "Zeggini.Pvalue\n");
    result.writeHeaderLine(fp);
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    result.clearValue();
    if (fitOK) {
      if (isBinaryOutcome()) {
        result.updateValue("Zeggini.Pvalue", logistic.GetPvalue());
      } else {
        result.updateValue("Zeggini.Pvalue", linear.GetPvalue());
      }
    }
    result.writeValueLine(fp);

    /* if (!fitOK) { */
    /*   fprintf(fp, "NA\n"); */
    /* } else { */
    /*   if (isBinaryOutcome()) { */
    /*     fprintf(fp, "%f\n", logistic.GetPvalue()); */
    /*   } else { */
    /*     fprintf(fp, "%f\n", linear.GetPvalue()); */
    /*   } */
    /* }; */
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
    this->modelName = "MadsonBrowning";
    result.addHeader("MB.Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();

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

    copyPhenotype(phenotype, &this->pheno);

    madsonBrowningCollapse(genotype, pheno, &collapsedGenotype);

    fitOK = logistic.FitNullModel(cov, pheno, 100);
    if (!fitOK) return -1;
    fitOK = logistic.TestCovariate(cov, pheno, collapsedGenotype);
    if (!fitOK) return -1;

    // record observed stat
    this->perm.init(logistic.GetStat()); // a chi-dist

    int failed = 0;
    while (this->perm.next()){
      permute(&this->pheno);
      madsonBrowningCollapse(genotype, pheno, &collapsedGenotype);
      fitOK = logistic.TestCovariate(collapsedGenotype, pheno);
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
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);

    if (isBinaryOutcome()){
      perm.writeHeaderTab(fp);
    } else {
      fp->write("Pvalue\n");
    };
    fp->write("\n");
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);

    if (isBinaryOutcome()) {
      perm.writeOutput(fp);
      fp->write("\n");
    } else {
      fp->write("NA\n");
    }
    /* if (isBinaryOutcome()) { */
    /*   if (fitOK){ */
    /*     perm.writeOutput(fp); */
    /*     fprintf(fp, "\n"); */
    /*   } else { */
    /*     fprintf(fp, "NA\tNA\tNA\tNA\tNA\tNA\n"); */
    /*   } */
    /* } else { */
    /*   fprintf(fp, "NA\n"); */
    /* } */
  };
private:
  Matrix collapsedGenotype;
  Vector pheno;
  LogisticRegressionScoreTest logistic;
  bool fitOK;
  int numVariant;
  Permutation perm;
}; // MadsonBrowningTest

// Danyu Lin's method, using 1/sqrt(p(1-p)) as weight
// where p is estimated from all samples
class FpTest: public ModelFitter{
public:
  FpTest() {
    this->modelName = "Fp";
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();

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
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    //fprintf(fp, "Fp.Pvalue\n");
    result.addHeader("Fp.Pvalue");
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (fitOK) {
      if (isBinaryOutcome()) {
        result.updateValue("Fp.Pvalue", logistic.GetPvalue());
      } else {
        result.updateValue("Fp.Pvalue", linear.GetPvalue());
      }
    }
    result.writeValueLine(fp);
    /* if (!fitOK) { */
    /*   fprintf(fp, "NA\n"); */
    /* } else { */
    /*   if (isBinaryOutcome()) { */
    /*     fprintf(fp, "%f\n", logistic.GetPvalue()); */
    /*   } else { */
    /*     fprintf(fp, "%f\n", linear.GetPvalue()); */
    /*   } */
    /* }; */
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
    this->result.addHeader("NumIncludeMarker");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();

    if (!isBinaryOutcome()) {
      fitOK = false;
      return -1;
    }
    if (covariate.cols != 0) { // rare cover does not take covariate
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
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.writeHeaderTab(fp);
    this->perm.writeHeaderLine(fp);

    /* if (isBinaryOutcome()){ */
    /*   fprintf(fp, "NumIncMarker\t"); */
    /*   this->perm.writeHeader(fp); */
    /*   fprintf(fp, "\n"); */
    /* } else */
    /*   fprintf(fp, "NA\n"); */
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (fitOK) {
      if (isBinaryOutcome()) {
        result.updateValue("NumIncludeMarker", (int)this->selected.size());
      }
    }
    result.writeValueTab(fp);
    this->perm.writeOutputLine(fp);


    /* if (isBinaryOutcome()) { */
    /*   if (fitOK){ */
    /*     fprintf(fp, "%zu\t", this->selected.size()); */
    /*     this->perm.writeOutput(fp); */
    /*     fprintf(fp, "\n"); */
    /*   } else { */
    /*     fprintf(fp, "NA\t"); */
    /*     this->perm.writeOutput(fp); */
    /*     fprintf(fp, "\n"); */
    /*   } */
    /* } else { */
    /*   fprintf(fp, "NA\n"); */
    /* } */
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

    double stat = -1;
    while ((int)selected.size() < genotype.rows) {
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
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate= dc->getCovariate();

    if (!isBinaryOutcome()) {
      fitOK = false;
      return -1;
    }
    if (covariate.cols != 0){
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
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    if (isBinaryOutcome()) { /// cmat only takes binary output
      this->perm.writeHeaderLine(fp);
    }
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (isBinaryOutcome()) {
      this->perm.writeOutputLine(fp);
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

    N_A = 0;
    N_U = 0;
    m_A = 0;
    m_U = 0;
    M_A = 0;
    M_U = 0;

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
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate= dc->getCovariate();
    Vector& weight = dc->getWeight();

    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    if (covariate.cols != 0) {
      fprintf(stderr, "Price's Variable Threshold method does not take covariate.\n");
    }

    rearrangeGenotypeByFrequency(genotype, &sortedGenotype, &this->freq);
    convertToReferenceAlleleCount(&sortedGenotype);
    collapseGenotype(&sortedGenotype);
    transpose(&sortedGenotype); // now each row is a collapsed genoype at certain frequency cutoff
    copyPhenotype(phenotype, &this->phenotype);

    this->zmax = -999.0;
    double z;
    if (isBinaryOutcome()) {
      for (int i = 0; i < sortedGenotype.rows; ++i) {
        z = calculateZForBinaryTrait(this->phenotype, this->sortedGenotype[i], weight);
        if ( z > this->zmax) {
          zmax = z;
          this->optimalFreq = freq[i];
        }
      }

      this->perm.init(zmax);

      // begin permutation
      while (this->perm.next()) {
        permute(&this->phenotype);
        double zp = -999.0;
        for (int j = 0; j < sortedGenotype.rows; ++j){
          z = calculateZForBinaryTrait(this->phenotype, this->sortedGenotype[j], weight);
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
        z = calculateZForContinuousTrait(this->phenotype, this->sortedGenotype[i], weight);
        if ( z > this->zmax){
          this->zmax = z;
          this->optimalFreq = freq[i];
        }
      }
      this->perm.init(this->zmax);

      // begin permutation
      while(this->perm.next()) {
        double zp = -999.0;
        permute(&this->phenotype);
        for (int j = 0; j < sortedGenotype.rows; ++j){
          z = calculateZForContinuousTrait(this->phenotype, this->sortedGenotype[j], weight);
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
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    fp->write("\tOptFreq\t");
    this->perm.writeHeader(fp);
    fp->write("\n");
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    fp->printf("\t%g\t", this->optimalFreq);
    this->perm.writeOutputLine(fp);
    // fprintf(fp, "\n");
  };
  void reset() {
    fitOK = true;
    this->perm.reset();
  };
private:
  /**
   * Convert genotype back to reference allele count
   * e.g. genotype 2 means homAlt/homAlt, so it has reference allele count 0
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
  double calculateZForBinaryTrait(Vector& y, Vector& x, Vector& weight){
    double ret = 0;
    int n = y.Length();

    if (weight.Length() == 0) {
      for (int i = 0; i < n; ++i) {
        if (y[i]) ret += x[i];
      }
    } else {
      for (int i = 0; i < n; ++i) {
        if (y[i]) ret += (x[i] * weight[i]);
      }
    }
    return ret;
  };
  double calculateZForContinuousTrait(Vector& y, Vector& x, Vector& weight){
    double ret = 0;
    int n = y.Length();
    if (weight.Length() == 0) {
      for (int i = 0; i < n; ++i) {
        ret += x[i] * y[i];
      }
    } else {
      for (int i = 0; i < n; ++i) {
        ret += x[i] * y[i] * weight[i];
      }
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
    result.addHeader("FreqThreshold");
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
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate= dc->getCovariate();

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
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    // this->result = siteInfo;
    // result.writeHeaderTab(fp);
    siteInfo.writeHeaderTab(fp);
    model[0].writeHeader(fp, result);
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    // char buf[1000];
    // siteInfo.writeValue(fp);
    for (size_t i = 0; i < freq.size(); i ++ ){
      /* sprintf(buf, "%s\t%f\t", prependString, freq[i]); */
      /* model[i].writeOutput(fp, buf); */
      result.updateValue("FreqThreshold", toString(freq[i]));
      //result.writeValueTab(fp);
      siteInfo.writeValueTab(fp);
      model[i].writeOutput(fp, result);
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
  // Result result;
}; // VariableThresholdCMCTest

#if 0
class VariableThresholdFreqTest: public ModelFitter{
public:
  // write result header
  void writeHeader(FileWriter* fp) {
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
  void writeOutput(FileWriter* fp) {
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
SkatTest(int nPerm, double alpha, double beta1, double beta2):perm(nPerm, alpha) {
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
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate= dc->getCovariate();
    Vector& weight = dc->getWeight();

    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
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
    this->res.Dimension(ynull.Length());
    for (int i = 0; i < this->ynull.Length(); ++i){
      this->res[i] = phenoVec[i] - this->ynull[i];
    }

    // get Pvalue
    skat.Fit(res, v, cov, genotype, weight);
    this->stat =  skat.GetQ();
    this->pValue = skat.GetPvalue();

    // permuation part
    this->perm.init(this->stat);

    double s;
    while (this->perm.next()) {
      permute(&res);
      s = skat.GetQFromNewResidual(res);
      this->perm.add(s);
    };
    fitOK = true;
    return 0;
  };
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    if (!usePermutation)
      fp->write("Q\tPvalue\n");
    else {
      fp->write("Q\tPvalue\t");
      this->perm.writeHeader(fp);
      fp->write("\n");
    }
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (!fitOK){
      fp->write("NA\tNA");
      if (usePermutation) {
        fp->write("\tNA\tNA\tNA\tNA\tNA\tNA");
      };
      fp->write("\n");
    } else {
      // binary outcome and quantative trait are similar output
      fp->printf("%g\t%g", this->skat.GetQ(), this->pValue);
      if (usePermutation) {
        fp->write("\t");
        this->perm.writeOutput(fp);
      }
      fp->write("\n");
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
  Vector res; //residual under the null
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
KbacTest(int nPerm, double alpha):nPerm(nPerm), alpha(alpha),
      xdatIn(NULL), ydatIn(NULL),mafIn(NULL), xcol(0), ylen(0), nn(0), qq(0) {
    this->modelName = "KBAC";
  };
  ~KbacTest() {
    if (this->xdatIn) delete[] this->xdatIn;
    if (this->ydatIn) delete[] this->ydatIn;
    if (this->mafIn) delete[] this->mafIn;
  };
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    fp->write("KBAC.Pvalue\n");
  };
  void reset() {
    // clear_kbac_test();
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate= dc->getCovariate();

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
    this->aa = this->alpha;
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
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (!fitOK){
      fp->write("NA\n");
    } else {
      if (isBinaryOutcome() ) {
        fp->printf("%f\n", this->pValue);
      } else {
        fp->printf("%f\n", this->pValue);
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
  double alpha;

  double* xdatIn;
  double* ydatIn;
  double* mafIn;
  int xcol;
  int ylen;
  int nn;
  int qq;
  double aa;
  double mafUpper;

  int twosided;
  bool fitOK;
  double pValue;
}; // KbacTest

//////////////////////////////////////////////////////////////////////
// Meta-analysis based methods
//////////////////////////////////////////////////////////////////////

// output files for meta-analysis
class MetaScoreTest: public ModelFitter{
public:
  MetaScoreTest(){
    this->modelName = "MetaScore";
  };
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate= dc->getCovariate();

    dc->countRawGenotype(0, &homRef, &het, &homAlt, &missing);
    // dc->getResult().writeValueLine(stderr);
    // fprintf(stderr, "%d\t%d\t%d\t%d\n", homRef, het, homAlt, missing);
    int nSample = (homRef + het + homAlt + missing);
    callRate = 1.0 - 1.0 * missing / nSample;
    if (homRef + het + homAlt == 0 ||
        (het < 0 || homRef < 0 || homAlt < 0)) {
      hweP = 0.0;
      af = 0.0;
    } else {
      hweP = SNPHWE( het, homRef, homAlt);
      af = 0.5 * (het + 2*homAlt) / (homRef + het + homAlt);
    }

    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }
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
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    if (g_SummaryHeader) {
      g_SummaryHeader->outputHeader(fp);
    }

    siteInfo.writeHeaderTab(fp);
    // fprintf(fp, "AF\tStat\tDirection\tPvalue\n");
    result.addHeader("AF");
    result.addHeader("INFORMATIVE_ALT_AC");
    result.addHeader("CALL_RATE");
    result.addHeader("HWE_PVALUE");
    result.addHeader("N_REF");
    result.addHeader("N_HET");
    result.addHeader("N_ALT");
    result.addHeader("U_STAT");
    result.addHeader("SQRT_V_STAT");
    result.addHeader("ALT_EFFSIZE");
    result.addHeader("PVALUE");
    result.writeHeaderLine(fp);
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {

    siteInfo.writeValueTab(fp);
    int informativeAC = het + 2* homAlt;

    result.clearValue();
    result.updateValue("AF", af);
    result.updateValue("INFORMATIVE_ALT_AC", informativeAC);
    result.updateValue("CALL_RATE", callRate);
    result.updateValue("HWE_PVALUE", hweP);
    result.updateValue("N_REF", homRef);
    result.updateValue("N_HET", het);
    result.updateValue("N_ALT", homAlt);
    if (fitOK) {
      if (!isBinaryOutcome()) {
        result.updateValue("U_STAT", linear.GetU()[0][0]);
        result.updateValue("SQRT_V_STAT", sqrt(linear.GetV()[0][0]));
        result.updateValue("ALT_EFFSIZE", linear.GetU()[0][0] / (linear.GetV()[0][0]));
        result.updateValue("PVALUE", linear.GetPvalue());
      } else {
        result.updateValue("U_STAT", logistic.GetU()[0][0]);
        result.updateValue("SQRT_V_STAT", sqrt(logistic.GetV()[0][0]));
        result.updateValue("ALT_EFFSIZE", logistic.GetU()[0][0] / (logistic.GetV()[0][0]));
        result.updateValue("PVALUE", logistic.GetPvalue());
      }
    }
    result.writeValueLine(fp);
  };
private:
  double af;
  // int nSample;
  Vector pheno;
  LinearRegressionScoreTest linear;
  LogisticRegressionScoreTest logistic;
  bool fitOK;
  Matrix cov;
  int homRef;
  int het;
  int homAlt;
  int missing;
  double hweP;
  double callRate;
  // Result result;
}; // MetaScoreTest

class MetaCovTest: public ModelFitter{
private:
  typedef std::vector<double> Genotype;
  struct Pos{
    std::string chrom;
    int pos;
  };
  struct Loci{
    Pos pos;
    Genotype geno;
  };

public:
  MetaCovTest(int windowSize){
    this->modelName = "MetaCov";
    this->numVariant = 0;
    this->nSample = -1;
    this->mleVarY = -1.;
    this->fout = NULL;
    this->windowSize = windowSize;
    result.addHeader("CHROM");
    result.addHeader("START_POS");
    result.addHeader("END_POS");
    result.addHeader("NUM_MARKER");
    result.addHeader("MARKER_POS");
    result.addHeader("COV");
  };
  ~MetaCovTest(){
    while(queue.size() > 0 ) {
      printCovariance(fout, queue);
      queue.pop_front();
    }
  }

  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Result& siteInfo = dc->getResult();

    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }

    if (genotype.rows == 0){
      fitOK = false;
      return -1;
    }
    if (nSample < 0) {
      nSample = genotype.rows;

      double s = 0;
      double s2 = 0;
      for (int i = 0; i < nSample; ++i) {
        s += phenotype[i][0];
        s2 += phenotype[i][0] * phenotype[i][0];
      }
      mleVarY = (s2 - s * s / nSample) / nSample;
      // fprintf(stderr, "MLE estimation of residual^2 = %g", mleVarY);
    } else {
      if (nSample != genotype.rows){
        fprintf(stderr, "Sample size changed at [ %s:%s ]", siteInfo["CHROM"].c_str(), siteInfo["POS"].c_str());
        fitOK = false;
        return -1;
      }
    }
    loci.pos.chrom = siteInfo["CHROM"];
    loci.pos.pos = atoi(siteInfo["POS"]);

    if ((siteInfo["REF"]).size() != 1 ||
        (siteInfo["ALT"]).size() != 1) { // not snp
      fitOK = false;
      return -1;
    };
    loci.geno.resize(nSample);
    for (int i = 0; i < nSample; ++i) {
      loci.geno[i] = genotype[i][0];
    }
    fitOK = true;
    return (fitOK ? 0 : -1);
  };

  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    if (g_SummaryHeader) {
      g_SummaryHeader->outputHeader(fp);
    }

    /* siteInfo.writeHeaderTab(fp); */
    // fprintf(fp, "AF\tStat\tDirection\tPvalue\n");
    result.writeHeaderLine(fp);
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    this->fout = fp;
    while (queue.size() && getWindowSize(queue, loci) > windowSize) {
      printCovariance(fout, queue);
      queue.pop_front();
    };
    if (fitOK) {
      queue.push_back(loci);
      ++numVariant;
    }
    // result.writeValueLine(fp);
  };

private:

  /**
   * @return \sum g1 * g2 - \sum(g1) * \sum(g2)/n
   */
  double getCovariance(const Genotype& g1, const Genotype& g2) {
    double sum_i = 0.0 ; // sum of genotype[,i]
    double sum_ij = 0.0 ; // sum of genotype[,i]*genotype[,j]
    double sum_j = 0.0 ; // sum of genotype[,j]
    int n = 0;
    for (size_t c = 0; c < g1.size(); ++c) { //iterator each people
      if (g1[c] < 0 || g2[c] < 0) continue;
      ++n;
      sum_i += g1[c];
      sum_ij += g1[c]*g2[c];
      sum_j += g2[c];
    };
    // fprintf(stderr, "n = %d sum_ij = %g sum_i = %g sum_j = %g \n", n, sum_ij, sum_i, sum_j);
    double cov_ij = (sum_ij - sum_i * sum_j / n) / n;
    // fprintf(stderr, "cov = %g var_i = %g var_j = %g n= %d\n", cov_ij, var_i, var_j, n);
    return cov_ij;
  };

  /**
   * @return max integer if different chromosome; or return difference between head and tail locus.
   */
  int getWindowSize(const std::deque<Loci>& loci, const Loci& newOne){
    if (loci.size() == 0) {
      return 0;
    }

    const Loci& head = loci.front();
    const Loci& tail = newOne;

    if (head.pos.chrom != tail.pos.chrom) {
      return INT_MAX;
    } else {
      return abs(tail.pos.pos - head.pos.pos);
    }
  };
  /**
   * @return 0
   * print the covariance for the front of loci to the rest of loci
   */
  int printCovariance(FileWriter* fp, const std::deque<Loci>& lociQueue){
    auto iter = lociQueue.begin();
    std::vector<int> position( lociQueue.size());
    std::vector<double> cov (lociQueue.size());
    int idx = 0;
    for (; iter != lociQueue.end(); ++iter){
      position[idx] = iter->pos.pos;
      cov[idx] = getCovariance(lociQueue.front().geno, iter->geno) * this-> mleVarY;
      idx ++;
    };
    /* fprintf(stderr, "%s\t%d\t%d\t", lociQueue.front().pos.chrom.c_str(), lociQueue.front().pos.pos, lociQueue.back().pos.pos);  */
    result.updateValue("CHROM", lociQueue.front().pos.chrom);
    result.updateValue("START_POS", lociQueue.front().pos.pos);
    result.updateValue("END_POS", lociQueue.back().pos.pos);
    /* fprintf(fp, "%d\t", idx); */
    result.updateValue("NUM_MARKER", idx);
    std::string s;
    for(int i = 0; i < idx; ++i) {
      if (i) s += ',';
      s+= toString(position[i]);
    }
    result.updateValue("MARKER_POS", s);

    s.clear();
    for(int i = 0; i < idx; ++i) {
      if (i) s+=',';
      s += floatToString(cov[i]);
    }
    result.updateValue("COV", s);
    result.writeValueLine(fp);
    return 0;
  };

private:
  std::deque< Loci> queue;
  int numVariant;
  int nSample;
  double mleVarY;
  FileWriter* fout;
  int windowSize;
  Loci loci;
  bool fitOK;
  // Result result;
}; // MetaCovTest


class DumpModel: public ModelFitter{
public:
  DumpModel(const char* prefix) {
    this->prefix = prefix;
    this->modelName = "DumpData";
  };
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) { // e.g. column headers.
    siteInfo.writeHeaderTab(fp);
    fp->write("FileName\n");

    this->header = siteInfo.joinHeader();
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate= dc->getCovariate();

    this->phenotype = phenotype;
    this->genotype = genotype;
    this->covariate = covariate;
    // copy covaraite column label
    for (int i = 0; i < covariate.cols; ++i){
      this->covariate.SetColumnLabel(i, covariate.GetColumnLabel(i));
    }
    // copy row labels
    this->rowLabel = dc->getRowLabel();
    return 0;
  };

  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo){
    std::string fn = this->prefix + "\t" + siteInfo.joinValue() + ".data";


    for (size_t i = 0; i < fn.size(); ++i) {
      if (fn[i] == '\t') fn[i] = '.';
    }

    siteInfo.writeValueTab(fp);
    fp->printf("%s\n", fn.c_str());

    // write header
    FILE* fDump= fopen(fn.c_str(), "wt");
    fprintf(fDump, "ID\t");
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
      if (strlen(covariate.GetColumnLabel(i)) == 0) {
        fprintf(fDump, "\tC%d", i);
      } else{
        fprintf(fDump, "\t%s", covariate.GetColumnLabel(i));
      }
    };
    fprintf(fDump, "\n");

    // write content
    for (int i = 0; i < phenotype.rows; ++i) {
      fprintf(fDump, "%s\t", rowLabel[i].c_str());
      // fputs(prependString, fDump);
      siteInfo.writeValue(fDump);
      for (int j = 0; j < phenotype.cols; ++j) {
        fprintf(fDump, "\t%f", phenotype[i][j]);
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
  std::vector<std::string> rowLabel;
}; // end DumpModel

//////////////////////////////////////////////////////////////////////
// Implementation of various collpasing methods
/**
 * @return Madson-Browning definition of alleleFrequency
 */
double getMarkerFrequency(Matrix& in, int col){
  int& numPeople = in.rows;
  double ac = 0; // NOTE: here genotype may be imputed, thus not integer
  int an = 0;
  for (int p = 0; p < numPeople; p++) {
    if ( (int) in[p][col] >= 0) {
      ac += in[p][col];
      an += 2;
    }
  }
  if ( an == 0 ) return 0.0;
  //double freq = 1.0 * (ac + 1) / (an + 1);
  double freq = ac / an;
  return freq;
};

double getMarkerFrequencyFromControl(Matrix& in, Vector& pheno, int col){
  int& numPeople = in.rows;
  double ac = 0; // NOTE: here genotype may be imputed, thus not integer
  int an = 0;
  for (int p = 0; p < numPeople; p++) {
    if (pheno[p] == 1) continue;
    if (in[p][col] >= 0) {
      ac += in[p][col];
      an += 2;
    }
  }
  // Refer:
  // 1. Madsen BE, Browning SR. A Groupwise Association Test for Rare Mutations Using a Weighted Sum Statistic. PLoS Genet. 2009;5(2):e1000384. Available at: http://dx.doi.org/10.1371/journal.pgen.1000384 [Accessed November 24, 2010].
  double freq = 1.0 * (ac + 1) / (an + 2);
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
    double weight = 1.0 / sqrt(freq * (1.0-freq)*genotype.rows);
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


#endif /* _MODELFITTER_H_ */
