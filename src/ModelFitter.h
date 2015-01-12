#ifndef _MODELFITTER_H_
#define _MODELFITTER_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <deque>

#include "libsrc/MathMatrix.h"

#include "base/ParRegion.h"
#include "regression/LogisticRegression.h"
#include "regression/LogisticRegressionScoreTest.h"
#include "regression/LogisticRegressionVT.h"
#include "regression/LinearRegression.h"
#include "regression/LinearRegressionScoreTest.h"
#include "regression/LinearRegressionVT.h"
#include "regression/Skat.h"
#include "regression/Table2by2.h"
#include "regression/kbac_interface.h"
#include "regression/FastLMM.h"
#include "regression/GrammarGamma.h"
#include "regression/MetaCov.h"
#include "regression/MatrixOperation.h"
#include "regression/MultivariateVT.h"
#include "regression/FamSkat.h"
#include "regression/FirthRegression.h"

#include "DataConsolidator.h"
#include "LinearAlgebra.h"
#include "ModelUtil.h"
#include "snp_hwe.c"
#include "Result.h"
#include "Summary.h"
#include "Permutation.h"

#if 0
// may decrease speed.
#ifdef _OPENMP
#include <omp.h>
#pragma message "Enable multithread using OpenMP"
#endif
#endif

extern SummaryHeader* g_SummaryHeader;

typedef void (*CollapsingFunction)(Matrix& in, const std::vector<int>& idx, Matrix* out, int index);

// various collapsing method
// they all take people by marker matrix
// and they won't take special care of missing genotypes
double getMarkerFrequency(Matrix& in, int col);
void getMarkerFrequency(Matrix& in, std::vector<double>* freq);
double getMarkerFrequencyFromControl(Matrix& in, Vector& pheno, int col);

void cmcCollapse(Matrix& in, Matrix* out);
void cmcCollapse(Matrix& in, const std::vector<int>& idx, Matrix* out, int index);

void zegginiCollapse(Matrix& in, Matrix* out);
void zegginiCollapse(Matrix& in, const std::vector<int>& idx, Matrix* out, int index);

void fpCollapse(Matrix& in, Matrix* out);

void madsonBrowningCollapse(Matrix& genotype, Vector& phenotype, Matrix* out);

void groupFrequency(const std::vector<double>& freq, std::map<double, std::vector<int> >* group);
void convertToMinorAlleleCount(Matrix& in, Matrix* g);
void convertToReferenceAlleleCount(Matrix& in, Matrix* g);

void makeVariableThreshodlGenotype(Matrix& in,
                                   const std::vector<double>& freqIn,
                                   Matrix* out,
                                   std::vector<double>* freqOut,
                                   void (*collapseFunc)(Matrix& , const std::vector<int>& , Matrix*, int)
                                   );

void appendHeritability(FileWriter* fp, const FastLMM& model);
void appendHeritability(FileWriter* fp, const GrammarGamma& model);

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
  // write footnotes in model output
  virtual void writeFootnote(FileWriter* fp) {};

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

class SingleVariantWaldTest: public ModelFitter{
 public:
  SingleVariantWaldTest(): fitOK(false){
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

class SingleVariantFirthTest: public ModelFitter{
 public:
  SingleVariantFirthTest(): fitOK(false){
    this->modelName = "SingleFirth";
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
      fitOK = false;
      return -1;
    }
    fitOK = firth.Fit(this->X, this->Y);
    // dumpToFile(X, "X");
    // dumpToFile(Y, "Y");
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
      result.clearValue();

      result.updateValue("Test", this->X.GetColumnLabel(i));
      if (fitOK) {
        result.updateValue("Beta", firth.GetCovEst()[i]);
        result.updateValue("SE", firth.GetCovB()[i][i]);
        result.updateValue("Pvalue", firth.GetAsyPvalue()[i]);
      }
    }
    result.writeValueLine(fp);
  }
 private:
  Matrix X; // 1 + cov + geno
  Vector Y; // phenotype
  FirthRegression firth;
  bool fitOK;
}; // SingleVariantFirthTest

class SingleVariantScoreTest: public ModelFitter{
 public:
  SingleVariantScoreTest(): af(-1), nSample(-1), fitOK(false),
                            needToFitNullModel(true){
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
      if (needToFitNullModel ||
          dc->isPhenotypeUpdated() ||
          dc->isCovariateUpdated()) {
        fitOK = linear.FitNullModel(cov, pheno);
        if (!fitOK) return -1;
        needToFitNullModel = false;
      }
      fitOK = linear.TestCovariate(cov, pheno, genotype);
    } else {
      if (needToFitNullModel ||
          dc->isPhenotypeUpdated() ||
          dc->isCovariateUpdated()) {
        fitOK = logistic.FitNullModel(cov, pheno, 100);
        if (!fitOK) return -1;
        needToFitNullModel = false;
      }
      fitOK = logistic.TestCovariate(cov, pheno, genotype);
    }
    return (fitOK ? 0 : -1);
  };

  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    /* if (g_SummaryHeader) { */
    /*   g_SummaryHeader->outputHeader(fp); */
    /* } */
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
  bool needToFitNullModel;
  Matrix cov;
}; // SingleVariantScoreTest

class SingleVariantFisherExactTest: public ModelFitter{
 public:
  SingleVariantFisherExactTest() {
    this->modelName = "FisherExact";
    caseAC = caseAN = ctrlAC = ctrlAN = 0;
    fitOK = false;
  };
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.addHeader("N00");
    result.addHeader("N01");
    result.addHeader("N10");
    result.addHeader("N11");
    result.addHeader("CtrlAF");
    result.addHeader("CaseAF");
    result.addHeader("PvalueTwoSide");
    result.addHeader("PvalueLess");
    result.addHeader("PvalueGreater");
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
      result.updateValue("N00", model.Get00());
      result.updateValue("N01", model.Get01());
      result.updateValue("N10", model.Get10());
      result.updateValue("N11", model.Get11());

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
      result.updateValue("PvalueTwoSide", model.getPExactTwoSided());
      result.updateValue("PvalueLess", model.getPExactOneSidedLess());
      result.updateValue("PvalueGreater", model.getPExactOneSidedGreater());
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
  SingleVariantFamilyScore():model(FastLMM::SCORE, FastLMM::MLE),
                             fitOK(false),
                             af(-1.), u(-1.), v(-1.), pvalue(-1.){
    this->modelName = "FamScore";
    result.addHeader("AF");
    result.addHeader("U.Stat");
    result.addHeader("V.Stat");
    result.addHeader("Pvalue");
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

    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }

    if (needToFitNullModel || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK = (0 == model.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(), *dc->getKinshipSForAuto()) ? true: false);
      if (!fitOK) return -1;
      needToFitNullModel = false;
    }

    fitOK = (0 == model.TestCovariate(cov, phenotype, genotype, *dc->getKinshipUForAuto(), *dc->getKinshipSForAuto()) ? true: false);
    af = model.GetAF(*dc->getKinshipUForAuto(), *dc->getKinshipSForAuto());
    u = model.GetUStat();
    v = model.GetVStat();
    pvalue = model.GetPvalue();
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
      if (isBinaryOutcome()) {

      } else {
        result.updateValue("AF", af);
        result.updateValue("U.Stat", u);
        result.updateValue("V.Stat", v);
        result.updateValue("Pvalue", pvalue);
      }
    }
    result.writeValueLine(fp);
  }
  void writeFootnote(FileWriter* fp) {
    appendHeritability(fp, model);
  }
 private:
  Matrix cov;
  FastLMM model;
  bool needToFitNullModel;
  bool fitOK;
  double af;
  double u;
  double v;
  double pvalue;
}; // end SingleVariantFamilyScore

class SingleVariantFamilyLRT: public ModelFitter{
 public:
  SingleVariantFamilyLRT():model(FastLMM::LRT, FastLMM::MLE),
                           fitOK(false),
                           af(-1.), nullLogLik(-1.), altLogLik(-1.), pvalue(-1.) {
    this->modelName = "FamLRT";
    result.addHeader("AF");
    result.addHeader("NullLogLik");
    result.addHeader("AltLogLik");
    result.addHeader("Pvalue");
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

    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }

    if (needToFitNullModel || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK = (0 == model.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(), *dc->getKinshipSForAuto()) ? true: false);
      if (!fitOK) return -1;
      needToFitNullModel = false;
    }

    fitOK = (0 == model.TestCovariate(cov, phenotype, genotype, *dc->getKinshipUForAuto(), *dc->getKinshipSForAuto()));
    af = model.GetAF(*dc->getKinshipUForAuto(), *dc->getKinshipSForAuto());
    nullLogLik = model.GetNullLogLikelihood();
    altLogLik = model.GetAltLogLikelihood();
    pvalue = model.GetPvalue();
    return (fitOK ? 0 : -1);
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.writeHeaderLine(fp);
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (fitOK) {
      if (isBinaryOutcome()) {
      } else {
        result.updateValue("AF", af);
        result.updateValue("NullLogLik", nullLogLik);
        result.updateValue("AltLogLik", altLogLik);
        result.updateValue("Pvalue", pvalue);
      }
    }
    result.writeValueLine(fp);
  }
  void writeFootnote(FileWriter* fp) {
    appendHeritability(fp, model);
  }

 private:
  Matrix cov;
  FastLMM model;
  bool needToFitNullModel;
  bool fitOK;
  double af;
  double nullLogLik;
  double altLogLik;
  double pvalue;
}; // end SingleVariantFamilyLRT

class SingleVariantFamilyGrammarGamma: public ModelFitter{
 public:
  SingleVariantFamilyGrammarGamma():model(),
                                    needToFitNullModel(true),
                                    fitOK(false),
                                    af (-1.),
                                    beta(-1.),
                                    betaVar(-1.),
                                    pvalue(-1.) {
    this->modelName = "FamGrammarGamma";
    result.addHeader("AF");
    result.addHeader("Beta");
    result.addHeader("BetaVar");
    result.addHeader("Pvalue");
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

    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }

    if (needToFitNullModel || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK = (0 == model.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(), *dc->getKinshipSForAuto()) ? true: false);
      if (!fitOK) return -1;
      needToFitNullModel = false;
    }

    fitOK = (0 == model.TestCovariate(cov, phenotype, genotype, *dc->getKinshipUForAuto(), *dc->getKinshipSForAuto()) ? true: false);
    af = model.GetAF(*dc->getKinshipUForAuto(), *dc->getKinshipSForAuto());
    beta = model.GetBeta();
    betaVar = model.GetBetaVar();
    pvalue = model.GetPvalue();
    return (fitOK ? 0 : -1);
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.writeHeaderLine(fp);
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (fitOK) {
      if (isBinaryOutcome()) {
      } else {
        result.updateValue("AF", af);
        result.updateValue("Beta", beta);
        result.updateValue("BetaVar", betaVar);
        result.updateValue("Pvalue", pvalue);
      }
    }
    result.writeValueLine(fp);
  }
  void writeFootnote(FileWriter* fp) {
    appendHeritability(fp, model);
  }

 private:
  Matrix cov;
  GrammarGamma model;
  bool needToFitNullModel;
  bool fitOK;
  double af;
  double beta;
  double betaVar;
  double pvalue;
}; // SingleVariantFamilyGrammarGamma

class CMCTest: public ModelFitter{
 public:
  CMCTest(): fitOK(false), numVariant(-1) {
    this->modelName = "CMC";
    result.addHeader("NonRefSite");
#if 0
    result.addHeader("U.Stat");
    result.addHeader("V.Stat");
    result.addHeader("Effect");
    result.addHeader("SE");
#endif
    result.addHeader("Pvalue");
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
#if 0
        result.updateValue("U.Stat", logistic.GetUStat());
        result.updateValue("V.Stat", logistic.GetVStat());
        result.updateValue("Effect", logistic.GetEffect());
        result.updateValue("SE"    , logistic.GetSE());
#endif
        result.updateValue("Pvalue", logistic.GetPvalue());
      } else {
#if 0
        result.updateValue("U.Stat", linear.GetUStat());
        result.updateValue("V.Stat", linear.GetVStat());
        result.updateValue("Effect", linear.GetEffect());
        result.updateValue("SE"    , linear.GetSE());
#endif
        result.updateValue("Pvalue", linear.GetPvalue());
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
  CMCWaldTest(): fitOK(false), numVariant(-1) {
    this->modelName = "CMCWald";
    result.addHeader("NonRefSite");
    result.addHeader("Beta");
    result.addHeader("SE");
    result.addHeader("Pvalue");
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
        result.updateValue("Pvalue", pval);
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

class ZegginiWaldTest: public ModelFitter{
 public:
  ZegginiWaldTest(): fitOK(false), numVariant(-1) {
    this->modelName = "ZegginiWald";
    result.addHeader("Beta");
    result.addHeader("SE");
    result.addHeader("Pvalue");
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

    zegginiCollapse(genotype, &collapsedGenotype);

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
        // result.updateValue("NonRefSite", this->totalNonRefSite());
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
        result.updateValue("Pvalue", pval);
      }
      result.writeValueLine(fp);
    }
  };
 private:
  Matrix collapsedGenotype;
  Matrix X;
  LogisticRegression logistic;
  LinearRegression linear;
  bool fitOK;
  int numVariant;
}; // ZegginiWaldTest

class CMCFisherExactTest: public ModelFitter{
 public:
  CMCFisherExactTest(): fitOK(false) {
    this->modelName = "CMCFisherExact";
    result.addHeader("N00");
    result.addHeader("N01");
    result.addHeader("N10");
    result.addHeader("N11");
    result.addHeader("PvalueTwoSide");
    result.addHeader("PvalueLess");
    result.addHeader("PvalueGreater");
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
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    result.clearValue();
    if (fitOK) {
      result.updateValue("N00", model.Get00());
      result.updateValue("N01", model.Get01());
      result.updateValue("N10", model.Get10());
      result.updateValue("N11", model.Get11());
      result.updateValue("PvalueTwoSide", model.getPExactTwoSided());
      result.updateValue("PvalueLess", model.getPExactOneSidedLess());
      result.updateValue("PvalueGreater", model.getPExactOneSidedGreater());

    }
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
  ZegginiTest(): fitOK(false), numVariant(-1) {
    this->modelName = "Zeggini";
    result.addHeader("Pvalue");
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
    //fprintf(fp, "Pvalue\n");
    result.writeHeaderLine(fp);
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    result.clearValue();
    if (fitOK) {
      if (isBinaryOutcome()) {
        result.updateValue("Pvalue", logistic.GetPvalue());
      } else {
        result.updateValue("Pvalue", linear.GetPvalue());
      }
    }
    result.writeValueLine(fp);
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
  MadsonBrowningTest(int nPerm, double alpha): fitOK(false),
                                               numVariant(-1),
                                               perm(nPerm, alpha) {
    this->modelName = "MadsonBrowning";
    result.addHeader("Pvalue");
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
    fitOK = false;
    numVariant = -1;
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
    //fprintf(fp, "Pvalue\n");
    result.addHeader("Pvalue");
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (fitOK) {
      if (isBinaryOutcome()) {
        result.updateValue("Pvalue", logistic.GetPvalue());
      } else {
        result.updateValue("Pvalue", linear.GetPvalue());
      }
    }
    result.writeValueLine(fp);
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
  RareCoverTest(int nPerm, double alpha): fitOK(false),
                                          numVariant(-1),
                                          stat(-1),
                                          perm(nPerm, alpha) {
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
    N_A = N_U = m_A = m_U = M_A = M_U = 0;
    fitOK = false;
    stat = -1.;
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
  VariableThresholdPrice(int nPerm, double alpha): fitOK(false),
                                                   zmax(-1.),
                                                   optimalFreq(-1),
                                                   perm(nPerm, alpha) {
    this->modelName = "VariableThresholdPrice";
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
    if (covariate.cols != 0) {
      fprintf(stderr, "Price's Variable Threshold method does not take covariate.\n");
    }

    // rearrangeGenotypeByFrequency(genotype, &sortedGenotype, &this->freq);
    this->freq.clear();
    makeVariableThreshodlGenotype(genotype, this->freq, &geno, &this->freq, cmcCollapse);
    convertToReferenceAlleleCount(geno, &sortedGenotype);
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
          }
          if (zp > this->zmax)  // early stop
            break;
        }
        this->perm.add(zp);
      }

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
          }
          if (zp > this->zmax) {
            break;
          }
        }
        this->perm.add(zp);
      }
    }

    fitOK = true;
    return 0;
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    fp->write("\tOptFreq\t");
    this->perm.writeHeader(fp);
    fp->write("\n");
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    fp->printf("\t%g\t", this->optimalFreq);
    this->perm.writeOutputLine(fp);
    // fprintf(fp, "\n");
  }
  void reset() {
    fitOK = true;
    this->perm.reset();
  }
 private:
  /**
   * @param g is people by marker matrix
   * we will collpase left to right, column by column
   * it is mimic the behavior of setting different frequency cutoff
   */
  void sequentialCollapseGenotype(Matrix* g){
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
  }
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
  }
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
  }
 private:
  Matrix geno;
  Matrix sortedGenotype;
  std::vector<double> freq;
  // std::map <double, std::vector<int> > freqGroup;

  bool fitOK;
  Vector phenotype;
  double zmax;
  double optimalFreq; // the frequency cutoff in unpermutated data which give smallest pvalue
  Permutation perm;
}; // VariableThresholdPrice


#if 0
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
    //rearrangeGenotypeByFrequency(genotype, &sortedGenotype, &this->freq);
    getMarkerFrequency(genotype, &freq);

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
#endif

class VTCMC: public ModelFitter{
 public:
  VTCMC(): fitOK(false), needToFitNullModel(true) {
    this->modelName = "VTCMC";
    result.addHeader("Freq");
    result.addHeader("U");
    result.addHeader("V");
    result.addHeader("OptimFreq");
    result.addHeader("Effect");
    result.addHeader("Pvalue");
  };
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate= dc->getCovariate();

    copy(phenotype, &this->pheno);
    // Vector& weight = dc->getWeight();

    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    // copyPhenotype(phenotype, &this->phenotype);
    int ret = fitNullModel(dc);
    if (ret) {
      return -1;
    }
    convertToMinorAlleleCount(genotype, &geno);
    // freq.clear();
    // for (int i = 0; i < genotype.rows; ++i){
    //   freq.push_back(getMarkerFrequency(geno, i));
    // }
    // groupFrequency(freq, &freqGroup);
    freq.clear();
    makeVariableThreshodlGenotype(geno, freq, &sortedGenotype, &freq, cmcCollapse);

    // rearrangeGenotypeByFrequency(genotype, &sortedGenotype, &this->freq);
    if (!isBinaryOutcome()) {
      ret = linear.TestCovariate(covariate, pheno, sortedGenotype);
    } else{
      ret = logistic.TestCovariate(covariate, pheno, sortedGenotype);
    }
    if (ret) {
      fitOK = false;
      return -1;
    }
    fitOK = true;
    return 0;
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.addHeader("OptFreq");
    result.writeHeaderLine(fp);
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    int index;
    double optimFreq;
    double effect;
    double pValue;
    if (!isBinaryOutcome()) {
      index = linear.GetIndexMax();
      optimFreq = freq[index];
      effect = linear.GetEffect(index);
      pValue = linear.GetPvalue();
      copyRowMatrix(linear.GetU(), &U);
      copyRowMatrix(linear.GetV(), &V);
    } else {
      index = logistic.GetIndexMax();
      optimFreq = freq[index];
      effect = logistic.GetEffect(index);
      pValue = logistic.GetPvalue();
      copyRowMatrix(logistic.GetU(), &U);
      copyRowMatrix(logistic.GetV(), &V);
    }
    result.updateValue("Freq", floatToString(freq));
    result.updateValue("U", floatToString(U));
    result.updateValue("V", floatToString(V));
    result.updateValue("OptimFreq", optimFreq);
    result.updateValue("Effect", effect);
    result.updateValue("Pvalue", pValue);

    result.writeValueLine(fp);
  };
  void reset() {
    fitOK = true;
  };

 private:
  int fitNullModel(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate= dc->getCovariate();
    if (this->needToFitNullModel || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &covariate);
      copyPhenotype(phenotype, &this->pheno);
      if (isBinaryOutcome()) {
        fitOK = logistic.FitNullModel(covariate, pheno);
      } else {
        fitOK = linear.FitNullModel(covariate, pheno);
      }
      if (!fitOK) return -1;
      needToFitNullModel = false;
    }
    return 0;
  }

  void copyRowMatrix(Matrix& m, std::vector<double>* out) {
    if (m.rows) return;
    int n = m.cols;
    out->resize(n);
    for (int i = 0; i < n; ++i) {
      (*out)[i] = m[0][i];
    }
  }
  Matrix geno;
  Matrix sortedGenotype;
  Vector pheno;
  std::vector<double> freq;
  std::vector<double> U;
  std::vector<double> V;
  LogisticRegressionVT logistic;
  LinearRegressionVT linear;
  bool fitOK;
  bool needToFitNullModel;
};

/**
 * Implementation of variable threshold from Liu's meta-analysis paper
 */
class AnalyticVT: public ModelFitter{
 public:
  typedef enum {
    UNRELATED = 0,
    RELATED = 1
  } Type;
 public:
  AnalyticVT(AnalyticVT::Type type):
      lmm(FastLMM::SCORE, FastLMM::MLE),
      fitOK(false),
      needToFitNullModel(true) {

    this->type = type;
    if (type == UNRELATED) {
      this->modelName = "AnalyticVT";
    } else {
      this->modelName = "FamAnalyticVT";
    }

    result.addHeader("MinMAF");
    result.addHeader("MaxMAF");
    result.addHeader("OptimMAF");
    result.addHeader("OptimNumVar");
    result.addHeader("U");
    result.addHeader("V");
    result.addHeader("Stat");
    result.addHeader("Pvalue");
  };
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate= dc->getCovariate();
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);

    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }

    if (isBinaryOutcome()) {
      // does not support binary outcomes
      return -1;
    }

    // obtain frequency, u and v per variant
    // const int nSample = genotype.rows;
    const int nVariant = genotype.cols;
    this->af.Dimension(nVariant);

    this->useFamilyModel = dc->hasKinship();
    if (this->useFamilyModel ^ (this->type == RELATED) ) {
      // model and data does not match
      return -1;
    }

    if (!this->useFamilyModel) {
      // calculate af
      for (int i = 0; i < nVariant; ++i) {
        af[i] = getMarkerFrequency(genotype, i);
      }

      // adjust covariates
      if (needToFitNullModel ||
          dc->isPhenotypeUpdated() ||
          dc->isCovariateUpdated()) {
        fitOK = linear.calculateResidualMatrix(cov, &residualMat);
        if (!fitOK) return -1;
        y.Product(residualMat, phenotype);
        needToFitNullModel = false;
      }
      x.Product(residualMat, genotype);
      centerMatrix(&x);

      // obtain sigma2
      sigma2 = getVariance(y, 0);

      // obtain U, V matrix
      xt.Transpose(x);
      u.Product(xt, y);
      v.Product(xt, x);
      v.Multiply(sigma2);
    } else {
      // family model
      if (needToFitNullModel ||
          dc->isPhenotypeUpdated() ||
          dc->isCovariateUpdated()) {
        fitOK = lmm.FitNullModel(cov, phenotype,
                                 *dc->getKinshipUForAuto(),
                                 *dc->getKinshipSForAuto()) == 0;
        if (!fitOK) return -1;
        needToFitNullModel = false;
      }

      for (int i = 0; i < nVariant; ++i) {
        af[i] = lmm.FastGetAF(*dc->getKinshipUForAuto(),
                              *dc->getKinshipSForAuto(),
                              genotype,
                              i);
        // fprintf(stderr, "af[%d] = %g\n", i, af[i]);
      }
      lmm.CalculateUandV(cov, phenotype, genotype,
                         *dc->getKinshipUForAuto(),
                         *dc->getKinshipSForAuto(),
                         &u, &v);
    }
    if (mvvt.compute(af, u, v)) {
      fitOK = false;
      return -1;
    }

    fitOK = true;
    return 0;
  }

  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.writeHeaderLine(fp);
  }

  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (fitOK) {
      result.updateValue("MinMAF", mvvt.getMinMAF());
      result.updateValue("MaxMAF", mvvt.getMaxMAF());
      result.updateValue("OptimMAF", mvvt.getOptimalMAF());
      result.updateValue("OptimNumVar", mvvt.getOptimalNumVar());
      result.updateValue("U", mvvt.getOptimalU());
      result.updateValue("V", mvvt.getOptimalV());
      result.updateValue("Stat", mvvt.getStat());
      result.updateValue("Pvalue", mvvt.getPvalue());
    }
    result.writeValueLine(fp);
  }
 private:
  Type type;
  Vector af;
  Matrix cov;
  Matrix x;
  Matrix xt;
  Matrix y;
  Matrix u;
  Matrix v;
  Matrix residualMat;
  double sigma2;
  LinearRegression linear;
  bool useFamilyModel;
  FastLMM lmm;
  MultivariateVT mvvt;
  bool fitOK;
  bool needToFitNullModel;
}; // AnalyticVT

class FamCMC: public ModelFitter{
 public:
  FamCMC(): needToFitNullModel(true),
            numVariant(0),
            u(-1.), v(-1.), af(-1.), effect(-1.), pvalue(-1.),
            lmm(FastLMM::SCORE, FastLMM::MLE),
            fitOK(false){
    this->modelName = "FamCMC";
    result.addHeader("NumSite");
    result.addHeader("AF");
    result.addHeader("U");
    result.addHeader("V");
    result.addHeader("Effect");
    result.addHeader("Pvalue");
  };
  // fitting model
  int fit(DataConsolidator* dc) {
    if (isBinaryOutcome()) {
      fitOK = false;
      return -1;
    }

    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();

    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }

    if (needToFitNullModel || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK = (0 == lmm.FitNullModel(cov, phenotype,
                                     *dc->getKinshipUForAuto(), *dc->getKinshipSForAuto()) ? true: false);
      if (!fitOK) return -1;
      needToFitNullModel = false;
    }

    cmcCollapse(genotype, &collapsedGenotype);

    dumpToFile(genotype, "genotype");
    dumpToFile(collapsedGenotype, "collapsedGenotype");

    fitOK = (0 == lmm.TestCovariate(cov, phenotype, collapsedGenotype,
                                    *dc->getKinshipUForAuto(),
                                    *dc->getKinshipSForAuto()));
    if (!fitOK) {
      return -1;
    }
    af = lmm.GetAF(*dc->getKinshipUForAuto(), *dc->getKinshipSForAuto(),
                   collapsedGenotype);
    u = lmm.GetUStat();
    v = lmm.GetVStat();
    if (v != 0) {
      effect = u / v;
    }
    pvalue = lmm.GetPvalue();
    return 0;
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.writeHeaderLine(fp);
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (fitOK && !isBinaryOutcome()) {
      result.updateValue("NumSite", numVariant);
      result.updateValue("AF", af);
      result.updateValue("U", u);
      result.updateValue("V", v);
      result.updateValue("Effect", effect);
      result.updateValue("Pvalue", pvalue);
    }
    result.writeValueLine(fp);
  }
  void writeFootnote(FileWriter* fp) {
    appendHeritability(fp, lmm);
  }
 private:
  Matrix cov;
  Matrix collapsedGenotype;
  bool needToFitNullModel;
  int numVariant;
  double u;
  double v;
  double af;
  double effect;
  double pvalue;
  FastLMM lmm;
  bool fitOK;
};

class FamZeggini: public ModelFitter{
 public:
  FamZeggini(): needToFitNullModel(true),
                numVariant(-1),
                u(-1.), v(-1.), af(-1.), effect(-1.), pvalue(-1.),
                lmm(FastLMM::SCORE, FastLMM::MLE),
                fitOK(false) {
    this->modelName = "FamZeggini";
    result.addHeader("NumSite");
    result.addHeader("MeanBurden");
    result.addHeader("U");
    result.addHeader("V");
    result.addHeader("Effect");
    result.addHeader("Pvalue");
  };
  // fitting model
  int fit(DataConsolidator* dc) {
    if (isBinaryOutcome()) {
      fitOK = false;
      return -1;
    }

    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();

    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }

    if (needToFitNullModel || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK = (0 == lmm.FitNullModel(cov, phenotype,
                                     *dc->getKinshipUForAuto(), *dc->getKinshipSForAuto()) ? true: false);
      if (!fitOK) return -1;
      needToFitNullModel = false;
    }

    zegginiCollapse(genotype, &collapsedGenotype);

    fitOK = (0 == lmm.TestCovariate(cov, phenotype, collapsedGenotype,
                                    *dc->getKinshipUForAuto(),
                                    *dc->getKinshipSForAuto()));
    if (!fitOK) {
      return -1;
    }
    af = lmm.GetAF(*dc->getKinshipUForAuto(), *dc->getKinshipSForAuto(),
                   collapsedGenotype);
    u = lmm.GetUStat();
    v = lmm.GetVStat();
    if (v != 0) {
      effect = u / v;
    }
    pvalue = lmm.GetPvalue();
    return 0;
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.writeHeaderLine(fp);
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (fitOK && !isBinaryOutcome()) {
      result.updateValue("NumSite", numVariant);
      result.updateValue("MeanBurden", af);
      result.updateValue("U", u);
      result.updateValue("V", v);
      result.updateValue("Effect", effect);
      result.updateValue("Pvalue", pvalue);
    }
    result.writeValueLine(fp);
  }
  void writeFootnote(FileWriter* fp) {
    appendHeritability(fp, lmm);
  }
 private:
  Matrix cov;
  Matrix collapsedGenotype;
  bool needToFitNullModel;
  int numVariant;
  double u;
  double v;
  double af;
  double effect;
  double pvalue;
  FastLMM lmm;
  bool fitOK;
};

class SkatTest: public ModelFitter{
 public:
  /* SkatTest(const std::vector<std::string>& param) { */
  SkatTest(int nPerm, double alpha, double beta1, double beta2):fitOK(false),
                                                                pValue(-1.),
                                                                stat(-1.),
                                                                perm(nPerm, alpha) {
    this->usePermutation = nPerm > 0;
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
    // fill it weight
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

  bool usePermutation;
  double stat;
  Permutation perm;
}; // SkatTest

class KBACTest: public ModelFitter{
 public:
  KBACTest(int nPerm, double alpha):nPerm(nPerm), alpha(alpha),
                                    xdatIn(NULL), ydatIn(NULL),mafIn(NULL), xcol(0), ylen(0), nn(0), qq(0),
                                    aa(alpha), mafUpper(1.), twosided(1), fitOK(false), pValue(-1.){
    this->modelName = "Kbac";
  };
  ~KBACTest() {
    if (this->xdatIn) delete[] this->xdatIn;
    if (this->ydatIn) delete[] this->ydatIn;
    if (this->mafIn) delete[] this->mafIn;
  };
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    fp->write("Pvalue\n");
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
      fp->printf("%f\n", this->pValue);
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
}; // KBACTest

class FamSkatTest: public ModelFitter{
 public:
  FamSkatTest(double beta1, double beta2):
      needToFitNullModel(true),
      fitOK(false),
      pValue(-1.),
      stat(-1.) {
    this->beta1 = beta1;
    this->beta2 = beta2;
    this->modelName = "FamSkat";
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

    if (isBinaryOutcome()) {
      fitOK = false;
      return -1;
    }
    const bool useFamilyModel = dc->hasKinship();
    if (!useFamilyModel) {
      fitOK = false;
      return -1;
    }

    // adjust covariates
    if (needToFitNullModel ||
        dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK = skat.FitNullModel(cov, phenotype,
                                *dc->getKinshipUForAuto(),
                                *dc->getKinshipSForAuto()) == 0;
      if (!fitOK) return -1;
      needToFitNullModel = false;
    }

    // get Pvalue
    fitOK = skat.TestCovariate(cov, phenotype, genotype, weight,
                               *dc->getKinshipUForAuto(),
                               *dc->getKinshipSForAuto()) == 0;
    if (!fitOK) {
      return -1;
    }

    this->stat =  skat.GetQ();
    this->pValue = skat.GetPvalue();

    fitOK = true;
    return 0;
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    fp->write("Q\tPvalue\n");
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (!fitOK){
      fp->write("NA\tNA");
      fp->write("\n");
    } else {
      // binary outcome and quantative trait are similar output
      fp->printf("%g\t%g", this->skat.GetQ(), this->pValue);
      fp->write("\n");
    }
  };
 private:
  bool needToFitNullModel;
  double beta1;
  double beta2;
  FamSkat skat;
  bool fitOK;
  double pValue;
  Matrix cov;
  double stat;
}; // FamSkatTest

//////////////////////////////////////////////////////////////////////
// Meta-analysis based methods
//////////////////////////////////////////////////////////////////////

// output files for meta-analysis
class MetaScoreTest: public ModelFitter{
 public:
  MetaScoreTest():
      linearFamScoreForAuto(FastLMM::SCORE, FastLMM::MLE),
      linearFamScoreForX(FastLMM::SCORE, FastLMM::MLE),
      needToFitNullModelForAuto(true),
      needToFitNullModelForX(true),
      hasHeritabilityForAuto(false),
      hasHeritabilityForX(false) {
    this->modelName = "MetaScore";
    af = afFromCase = afFromControl = -1.;
    fitOK = false;
    homRef = het = homAlt = missing = -1;
    hweP = hwePvalueFromCase = hwePvalueFromControl = -1.;
    callRate = -1.;
    useFamilyModel = false;
    isHemiRegion = false;
  }
  // fitting model
  virtual int fit(DataConsolidator* dc) {
    Matrix& genotype  = dc->getGenotype();
    return this->fitWithGivenGenotype(genotype, dc);
  }
  int fitWithGivenGenotype(Matrix& genotype, DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& covariate = dc->getCovariate();

    // check column name for hemi region
    this->isHemiRegion = dc->isHemiRegion(0);

    dc->countRawGenotype(0, &homRef, &het, &homAlt, &missing);
    // dc->getResult().writeValueLine(stderr);
    // fprintf(stderr, "%d\t%d\t%d\t%d\n", homRef, het, homAlt, missing);
    int nSample = (homRef + het + homAlt + missing);
    if (nSample){
      callRate = 1.0 - 1.0 * missing / nSample;
    } else {
      callRate = 0.0;
    }
    // use all samples to get af
    if (nSample){
      af = 0.5 * (het + 2. * homAlt) / (homRef + het + homAlt);
    } else {
      af = -1.0;
    }
    // fprintf(stderr, "af = %g\n", af);
    // use female info to get HWE
    if (isHemiRegion) {
      dc->countRawGenotypeFromFemale(0, &homRef, &het, &homAlt, &missing);
    }
    if (homRef + het + homAlt == 0 ||
        (het < 0 || homRef < 0 || homAlt < 0)) {
      hweP = 0.0;
      hwePvalueFromCase = 0.0;
      hwePvalueFromControl = 0.0;
      // af = 0.0;
      // afFromCase = 0.0;
      // afFromControl = 0.0;
    } else {
      hweP = SNPHWE( het, homRef, homAlt);
      // af = 0.5 * (het + 2*homAlt) / (homRef + het + homAlt);
    }

    // handle binary cases
    if (isBinaryOutcome()) {
      int homRefCase, hetCase, homAltCase, missingCase;
      int homRefCtrl, hetCtrl, homAltCtrl, missingCtrl;

      if (!isHemiRegion) {
        dc->countRawGenotypeFromCase(0, &homRefCase, &hetCase, &homAltCase, &missingCase);
        dc->countRawGenotypeFromControl(0, &homRefCtrl, &hetCtrl, &homAltCtrl, &missingCtrl);
      } else {
        dc->countRawGenotypeFromFemaleCase(0, &homRefCase, &hetCase, &homAltCase, &missingCase);
        dc->countRawGenotypeFromFemaleControl(0, &homRefCtrl, &hetCtrl, &homAltCtrl, &missingCtrl);
      }

      if (homRefCase + hetCase + homAltCase == 0 ||
          (hetCase < 0 || homRefCase < 0 || homAltCase < 0)) {
        hwePvalueFromCase = 0.0;
        afFromCase = 0.0;
      } else {
        hwePvalueFromCase = SNPHWE(hetCase, homRefCase, homAltCase);
        afFromCase = 0.5 * (hetCase + 2*homAltCase) / (homRefCase + hetCase + homAltCase);
      }

      if (homRefCtrl + hetCtrl + homAltCtrl == 0 ||
          (hetCtrl < 0 || homRefCtrl < 0 || homAltCtrl < 0)) {
        hwePvalueFromControl = 0.0;
        afFromControl = 0.0;
      } else {
        hwePvalueFromControl = SNPHWE(hetCtrl, homRefCtrl, homAltCtrl);
        afFromControl = 0.5 * (hetCtrl + 2*homAltCtrl) / (homRefCtrl + hetCtrl + homAltCtrl);
      }
    }

    // sanity check, this should not happen
    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }

    // skip monomorphic sites
    if (isMonomorphicMarker(genotype, 0)) {
      fitOK = false;
      return -1;
    }

    // perform assocation tests
    this->useFamilyModel = dc->hasKinship();
    if (this->useFamilyModel) {
      if (!isBinaryOutcome()) { // quant trait
        if (!isHemiRegion) {
          if (needToFitNullModelForAuto ||
              dc->isPhenotypeUpdated() ||
              dc->isCovariateUpdated()) {
            copyCovariateAndIntercept(genotype.rows, covariate, &cov);
            fitOK = (0 == linearFamScoreForAuto.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(), *dc->getKinshipSForAuto()) ? true: false);
            if (!fitOK) return -1;
            hasHeritabilityForAuto = true;
            needToFitNullModelForAuto = false;
          }
        } else { // hemi region
          if (!dc->hasKinshipForX()) {
            fitOK = false;
            return -1;
          }
          if (needToFitNullModelForX ||
              dc->isPhenotypeUpdated() ||
              dc->isCovariateUpdated()) {
            copyCovariateAndIntercept(genotype.rows, covariate, &cov);
            fitOK = (0 == linearFamScoreForX.FitNullModel(cov, phenotype, *dc->getKinshipUForX(), *dc->getKinshipSForX()) ? true: false);
            if (!fitOK) return -1;
            hasHeritabilityForX = true;
            needToFitNullModelForX = false;
          }
        }

        if (!isHemiRegion) {
          fitOK = (0 == linearFamScoreForAuto.TestCovariate(cov, phenotype, genotype, *dc->getKinshipUForAuto(), *dc->getKinshipSForAuto()) ? true: false);
          this->af = linearFamScoreForAuto.FastGetAF(*dc->getKinshipUForAuto(), *dc->getKinshipSForAuto(), dc->getGenotype());
        } else {
          fitOK = (0 == linearFamScoreForX.TestCovariate(cov, phenotype, genotype, *dc->getKinshipUForX(), *dc->getKinshipSForX()) ? true: false);
          this->af = linearFamScoreForX.FastGetAF(*dc->getKinshipUForX(), *dc->getKinshipSForX(), dc->getGenotype());
        }
      } else {
        /* if (needToFitNullModel || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) { */
        /*   copyCovariateAndIntercept(genotype.rows, covariate, &cov); */
        /*   copyPhenotype(phenotype, &this->pheno); */
        /*   fitOK = logistic.FitNullModel(cov, pheno, 100); */
        /*   if (!fitOK) return -1; */
        /*   needToFitNullModel = false; */
        /* } */
        /* fitOK = logistic.TestCovariate(cov, pheno, genotype); */
        // return -1;
        fprintf(stderr, "Cannot support binary related individuals\n");
        exit(1);
      }
    } else { // unrelated
      if (!isBinaryOutcome()) { // continuous/quantative trait
        if (this->needToFitNullModelForAuto || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) {
          copyCovariateAndIntercept(genotype.rows, covariate, &cov);
          copyPhenotype(phenotype, &this->pheno);
          fitOK = linear.FitNullModel(cov, pheno);
          if (!fitOK) return -1;
          needToFitNullModelForAuto = false;
        }
        fitOK = linear.TestCovariate(cov, pheno, genotype);
        if (!fitOK) return -1;
      } else { // binary trait
        // fit null model
        if (this->needToFitNullModelForAuto || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) {
          copyCovariateAndIntercept(genotype.rows, covariate, &cov);
          copyPhenotype(phenotype, &this->pheno);
          fitOK = logistic.FitNullModel(cov, pheno, 100);
          if (!fitOK) return -1;
          needToFitNullModelForAuto = false;
        }
        fitOK = logistic.TestCovariate(cov, pheno, genotype);
        if (!fitOK) return -1;
        // fit alternative model
        if (covariate.cols) {
          copyGenotypeWithCovariateAndIntercept(genotype, covariate, &this->X);
        } else {
          copyGenotypeWithIntercept(genotype, &this->X);
        }
        Vector& b_null = logistic.GetNullCovEst();
        Vector b(b_null.Length() + 1);
        for (int i = 0; i < b_null.Length(); ++i) {
          b[i] = b_null[i];
        }
        b[b_null.Length()] = 0.0;

        logisticAlt.SetInitialCovEst(b);
        fitOK = logisticAlt.FitLogisticModel(this->X, this->pheno, 100);
        if (!fitOK) return -1;
      }
    }
    return (fitOK ? 0 : -1);
  }

  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    if (g_SummaryHeader) {
      g_SummaryHeader->outputHeader(fp);
    }

    siteInfo.writeHeaderTab(fp);
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
    if (af > 0) {
      if (!isBinaryOutcome()) {
        result.updateValue("AF", af);
      } else {
        static char afString[128];
        snprintf(afString, 128, "%g:%g:%g", af, afFromCase, afFromControl);
        result.updateValue("AF", afString);
      }
    }

    result.updateValue("INFORMATIVE_ALT_AC", informativeAC);
    result.updateValue("CALL_RATE", callRate);
    if (!isBinaryOutcome()) {
      result.updateValue("HWE_PVALUE", hweP);
    } else {
      static char hwePString[128];
      snprintf(hwePString, 128, "%g:%g:%g", hweP, hwePvalueFromCase, hwePvalueFromControl);
      result.updateValue("HWE_PVALUE", hwePString);
    }
    result.updateValue("N_REF", homRef);
    result.updateValue("N_HET", het);
    result.updateValue("N_ALT", homAlt);
    if (fitOK) {
      if (this->useFamilyModel) {
        if (!isBinaryOutcome()) {
          if (!isHemiRegion) {
            result.updateValue("U_STAT", linearFamScoreForAuto.GetUStat());
            result.updateValue("SQRT_V_STAT", sqrt(linearFamScoreForAuto.GetVStat()));
            result.updateValue("ALT_EFFSIZE", linearFamScoreForAuto.GetVStat() != 0.0 ?
                               linearFamScoreForAuto.GetUStat() / linearFamScoreForAuto.GetVStat() : 0.0);
            result.updateValue("PVALUE", linearFamScoreForAuto.GetPvalue());
          } else {
            result.updateValue("U_STAT", linearFamScoreForX.GetUStat());
            result.updateValue("SQRT_V_STAT", sqrt(linearFamScoreForX.GetVStat()));
            result.updateValue("ALT_EFFSIZE", linearFamScoreForX.GetVStat() != 0.0 ?
                               linearFamScoreForX.GetUStat() / linearFamScoreForX.GetVStat() : 0.0);
            result.updateValue("PVALUE", linearFamScoreForX.GetPvalue());
          }
        } else {
          /* result.updateValue("U_STAT", logistic.GetU()[0][0]); */
          /* result.updateValue("SQRT_V_STAT", sqrt(logistic.GetV()[0][0])); */
          /* result.updateValue("ALT_EFFSIZE", logistic.GetU()[0][0] / (logistic.GetV()[0][0])); */
          /* result.updateValue("PVALUE", logistic.GetPvalue()); */
        }
      } else { // unrelated individual
        if (!isBinaryOutcome()) {
          const double sigma2 = linear.GetSigma2();
          if (sigma2 != 0.0) {
            // U => U / sigma^2
            // V => V / sigma^2 / sigma^2 (NOTE: V includes the sigma^2 term)
            result.updateValue("U_STAT", linear.GetU()[0][0] / sigma2);
            result.updateValue("SQRT_V_STAT", sqrt(linear.GetV()[0][0]) / sigma2 );
            result.updateValue("ALT_EFFSIZE", linear.GetV()[0][0] != 0.0 ? linear.GetBeta()[0][0] : 0.0);
            result.updateValue("PVALUE", linear.GetPvalue());
          }
        } else {
          result.updateValue("U_STAT", logistic.GetU()[0][0]);
          result.updateValue("SQRT_V_STAT", sqrt(logistic.GetV()[0][0]));
          result.updateValue("ALT_EFFSIZE", logisticAlt.GetCovEst()[1]); // first column is intercept, so 1 means second column which is beta
          result.updateValue("PVALUE", logistic.GetPvalue());
        }
      }
    }
    result.writeValueLine(fp);
  }

  void writeFootnote(FileWriter* fp) {
    if (this->hasHeritabilityForAuto) {
      fp->writeLine("#Region\tautosomal\n");
      appendHeritability(fp, linearFamScoreForAuto);
    }
    if (this->hasHeritabilityForX) {
      fp->writeLine("#Region\tX\n");
      appendHeritability(fp, linearFamScoreForX);
    }
  }

 private:
  double af;
  double afFromCase;
  double afFromControl;
  Vector pheno;

  // QT linear model for unrelated
  LinearRegressionScoreTest linear;
  // Binary logistic model for unrelated
  LogisticRegressionScoreTest logistic;
  LogisticRegression logisticAlt;
  // QT linear model for related
  FastLMM linearFamScoreForAuto;
  FastLMM linearFamScoreForX;

  bool fitOK;
  Matrix cov;
  Matrix X; // intercept, cov(optional) and genotype

  int homRef;
  int het;
  int homAlt;
  int missing;
  double hweP;
  double hwePvalueFromCase;
  double hwePvalueFromControl;
  double callRate;

  bool needToFitNullModelForAuto;
  bool needToFitNullModelForX;
  bool useFamilyModel;
  bool hasHeritabilityForAuto;
  bool hasHeritabilityForX;

  bool isHemiRegion; // is the variant tested in hemi region?
}; // MetaScoreTest

class MetaDominantTest: public MetaScoreTest{
 public:
  MetaDominantTest(): MetaScoreTest() {
    this->modelName = "MetaDominant";
  }
  virtual int fit(DataConsolidator* dc) {
    dc->codeGenotypeForDominantModel(&geno);
    return fitWithGivenGenotype(geno, dc);
  }
 private:
  Matrix geno;
};

class MetaRecessiveTest: public MetaScoreTest{
 public:
  MetaRecessiveTest(): MetaScoreTest() {
    this->modelName = "MetaRecessive";
  }
  virtual int fit(DataConsolidator* dc) {
    dc->codeGenotypeForRecessiveModel(&geno);
    return fitWithGivenGenotype(geno, dc);
  }
 private:
  Matrix geno;
};

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
  MetaCovTest(int windowSize): fitOK(false),
                               useFamilyModel(false),
                               isHemiRegion(false) {
    this->modelName = "MetaCov";
    this->numVariant = 0;
    this->nSample = -1;
    // this->mleVarY = -1.;
    this->fout = NULL;
    this->windowSize = windowSize;
    this->needToFitNullModelForAuto = true;
    this->needToFitNullModelForX = true;
    result.addHeader("CHROM");
    result.addHeader("START_POS");
    result.addHeader("END_POS");
    result.addHeader("NUM_MARKER");
    result.addHeader("MARKER_POS");
    result.addHeader("COV");
    /* result.addHeader("COV_XZ"); */
    /* result.addHeader("COV_ZZ"); */
  };
  ~MetaCovTest(){
    while(queue.size() > 0 ) {
      if (isBinaryOutcome()) {
        printCovarianceForBinaryTrait(fout, queue);
      } else {
        printCovariance(fout, queue);
      }
      queue.pop_front();
    }
  }

  // fitting model
  virtual int fit(DataConsolidator* dc) {
    Matrix& genotype  = dc->getGenotype();
    return this->fitWithGivenGenotype(genotype, dc);
  }
  int fitWithGivenGenotype(Matrix& genotype, DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& covariate = dc->getCovariate();
    Result& siteInfo = dc->getResult();
    this->isHemiRegion = dc->isHemiRegion(0);

    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }
    if (genotype.rows == 0){
      fitOK = false;
      return -1;
    }
    this->useFamilyModel = dc->hasKinship();
    if (nSample < 0) { // unitialized
      // calculate variance of y
      nSample = genotype.rows;
      weight.Dimension(nSample);
    }
    if (nSample != genotype.rows){
      fprintf(stderr, "Sample size changed at [ %s:%s ]", siteInfo["CHROM"].c_str(), siteInfo["POS"].c_str());
      fitOK = false;
      return -1;
    }

    // set weight
    if (!isBinaryOutcome()) { // qtl
      if (useFamilyModel) {
        if (!isHemiRegion) { // autosomal
          // fit null model
          if (this->needToFitNullModelForAuto || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) {
            copyCovariateAndIntercept(genotype.rows, covariate, &cov);
            fitOK = (0 == metaCovForAuto.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(), *dc->getKinshipSForAuto()));
            if (!fitOK) {
              return -1;
            }
            needToFitNullModelForAuto = false;
            // get weight
            metaCovForAuto.GetWeight(&this->weight);
          }
        } else { // hemi region
          if (!dc->hasKinshipForX()) {
            fitOK = false;
            return -1;
          }
          if (this->needToFitNullModelForX || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) {
            copyCovariateAndIntercept(genotype.rows, covariate, &cov);
            fitOK = (0 == metaCovForX.FitNullModel(cov, phenotype, *dc->getKinshipUForX(), *dc->getKinshipSForX()));
            if (!fitOK) {
              return -1;
            }
            needToFitNullModelForX = false;
            // get weight
            metaCovForX.GetWeight(&this->weight);
          }
        }
      } else { // not family model
        double s = 0;
        double s2 = 0;
        for (int i = 0; i < nSample; ++i) {
          s += phenotype[i][0];
          s2 += phenotype[i][0] * phenotype[i][0];
        }
        double sigma2 = (s2 - s * s / nSample) / nSample;
        // fprintf(stderr, "sigma2 = %g\n", sigma2);
        if (sigma2 != 0) {
          for (int i = 0; i < nSample ; ++i) {
            weight[i]  = 1.0 / sigma2;  // mleVarY
          }
        } else{
          fprintf(stderr, "sigma2 = 0.0 for the phenotype!\n");
          for (int i = 0; i < nSample ; ++i) {
            weight[i]  = 1.0;
          }
        }
      }
      // fprintf(stderr, "MLE estimation of residual^2 = %g", mleVarY);
    } else { // binary case
      if (useFamilyModel) {
        fprintf(stderr, "Binary trait meta covariance are not supported!\n");
        exit(1); // not supported yet
      } else {
        // fit null model, don't distinguish chrom X in this case
        if (this->needToFitNullModelForAuto || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) {
          copyCovariateAndIntercept(genotype.rows, covariate, &cov);
          copyPhenotype(phenotype, &this->pheno);
          fitOK = logistic.FitLogisticModel(cov, pheno, 100);
          if (!fitOK) return -1;
          needToFitNullModelForAuto = false;

          // skip store Z, as Z = this->cov
          // store V in weight
          for (int i = 0; i < nSample; ++i) {
            const double y = logistic.GetPredicted()[i];
            weight[i] = y * (1.0 - y);
          }
          // store Z^T * V * Z
          this->ZVZ.Dimension(cov.cols, cov.cols);
          for(int i = 0; i < cov.cols; ++i) {
            for (int j = 0; j < cov.cols; ++j) {
              double s = 0.0;
              for (int k = 0; k < cov.rows; ++k) {
                s += cov[k][i] * weight[k] * cov[k][j];
              }
              this->ZVZ[i][j] = s;
            }
          }
        }
      }
    }
    loci.pos.chrom = siteInfo["CHROM"];
    loci.pos.pos = atoi(siteInfo["POS"]);

    if ((siteInfo["REF"]).size() != 1 ||
        (siteInfo["ALT"]).size() != 1) { // not snp
      fitOK = false;
      return -1;
    };

    // assign loci.geno, and
    // check if this is a monomorphic site, if so, just skip it.
    loci.geno.resize(nSample);
    double g0 = genotype[0][0];
    bool isVarint = false;
    for (int i = 0; i < nSample; ++i) {
      loci.geno[i] = genotype[i][0];
      if (loci.geno[i] != g0)
        isVarint = true;
    }
    if (!isVarint) {
      fitOK = false;
      return -1;
    }

    if (!isBinaryOutcome()) {
      if (useFamilyModel) {
        if (!isHemiRegion) {
          metaCovForAuto.TransformCentered(&loci.geno,
                                           *dc->getKinshipUForAuto(),
                                           *dc->getKinshipSForAuto());
        } else {
          metaCovForX.TransformCentered(&loci.geno,
                                        *dc->getKinshipUForX(),
                                        *dc->getKinshipSForX());
        }
      } else { // unrelated individuals
        // center genotype
        double s = 0.0;
        for (int i = 0; i < nSample; ++i) {
          s += loci.geno[i];
        }
        s /= nSample;
        for (int i = 0; i < nSample; ++i) {
          loci.geno[i] -= s;
        }
      }
    } else {
      /// for binary, no need to center. and not support family.
    }
    fitOK = true;
    return (fitOK ? 0 : -1);
  }; // fitWithGivenGenotype

  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    if (g_SummaryHeader) {
      g_SummaryHeader->outputHeader(fp);
    }

    result.writeHeaderLine(fp);
  };
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    this->fout = fp;
    while (queue.size() && getWindowSize(queue, loci) > windowSize) {
      if (isBinaryOutcome()) {
        printCovarianceForBinaryTrait(fout, queue);
      } else {
        printCovariance(fout, queue);
      }
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
   * NOTE: we already centered g1, g2, so the above reduced to
   *        (U %*% (g1 - \bar(g1)*))' %*% weight %*% (U %*% (g2 - \bar(g2)))
   */
  double getCovariance(const Genotype& g1, const Genotype& g2) {
    double sum_ij = 0.0 ; // sum of genotype[,i]*genotype[,j]
    int n = g1.size();
    if (n == 0) return 0.0;
    for (int c = 0; c < n; ++c) { //iterator each people
      // fprintf(stderr, "weight[%d] = %g\n", (int)c, weight[int(c)]);
      sum_ij += g1[c]*g2[c] / this->weight[(int)c];
    }
    double cov_ij = sum_ij / n;
    
    return cov_ij;
  };

  /**
   * @return g1^T * weight * g2
   */
  double getCovarianceForBinaryTrait(const Genotype& g1, const Genotype& g2) {
    double sum = 0.0;
    const int n = g1.size();
    if (n == 0) return 0.0;
    for (int i = 0; i < n; ++i){
      sum += g1[i] * this->weight[(int)i] * g2[i];
    }
    return sum;
  };
  /**
   * @return g^T * weight * m[,col]
   */
  double getCovarianceForBinaryTrait(const Genotype& g, Matrix& m, int col) {
    double sum = 0.0;
    const int n = g.size();
    if (n == 0) return 0.0;
    for (int i = 0; i < n; ++i){
      sum += g[i] * this->weight[(int)i] * m[i][col];
    }
    return sum;
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
    std::deque<Loci>::const_iterator iter = lociQueue.begin();
    std::vector<int> position(lociQueue.size());
    std::vector<double> cov (lociQueue.size());
    int idx = 0;
    for (; iter != lociQueue.end(); ++iter){
      position[idx] = iter->pos.pos;
      // cov[idx] = getCovariance(lociQueue.front().geno, iter->geno) * this-> mleVarY;
      cov[idx] = getCovariance(lociQueue.front().geno, iter->geno);
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
  /**
   * @return 0
   * print the covariance for the front of loci to the rest of loci
   */
  int printCovarianceForBinaryTrait(FileWriter* fp, const std::deque<Loci>& lociQueue){
    std::deque<Loci>::const_iterator iter = lociQueue.begin();
    std::vector<int> position(lociQueue.size());
    std::vector<double> covXX (lociQueue.size());
    std::vector<double> covXZ (this->cov.cols);
    const size_t numMarker = lociQueue.size();

    for (int idx = 0; iter != lociQueue.end(); ++iter, ++idx){
      position[idx] = iter->pos.pos;
      covXX[idx] = getCovarianceForBinaryTrait(lociQueue.front().geno, iter->geno);
    }
    for (int i = 0; i < this->cov.cols; ++i) {
      covXZ[i] = getCovarianceForBinaryTrait(lociQueue.front().geno, this->cov, i);
    }
    /* fprintf(stderr, "%s\t%d\t%d\t", lociQueue.front().pos.chrom.c_str(), lociQueue.front().pos.pos, lociQueue.back().pos.pos);  */
    result.updateValue("CHROM", lociQueue.front().pos.chrom);
    result.updateValue("START_POS", lociQueue.front().pos.pos);
    result.updateValue("END_POS", lociQueue.back().pos.pos);
    /* fprintf(fp, "%d\t", idx); */
    result.updateValue("NUM_MARKER", (int)numMarker);
    std::string s;
    for(size_t i = 0; i < numMarker; ++i) {
      if (i) s += ',';
      s+= toString(position[i]);
    }
    result.updateValue("MARKER_POS", s);

    s.clear();
    for(size_t i = 0; i < numMarker; ++i) {
      if (i) s+=',';
      s += floatToString(covXX[i]);
    }
    s += ':';
    for(int i = 0; i < this->cov.cols; ++i) {
      if (i) s+=',';
      s += floatToString(covXZ[i]);
    }
    s += ':';
    for (int i = 0; i < this->ZVZ.rows; ++i) {
      for (int j = 0; j <= i ; ++j) {
        if (i || j) s+=',';
        s += floatToString(this->ZVZ[i][j]);
      }
    }
    result.updateValue("COV", s);
    result.writeValueLine(fp);
    return 0;
  };

 private:
  std::deque< Loci> queue;
  int numVariant;
  int nSample;
  // Vector mleVarY;  // variance term of Y_i for i = 1,..., N th sample
  FileWriter* fout;
  int windowSize;
  Loci loci;
  bool fitOK;
  // Result result;
  MetaCov metaCovForAuto;
  MetaCov metaCovForX;
  bool needToFitNullModelForAuto;
  bool needToFitNullModelForX;
  bool useFamilyModel;
  Vector weight; // per individual weight
  LogisticRegression logistic;
  Matrix cov;
  Matrix ZVZ;
  Vector pheno;
  bool isHemiRegion; // is the variant tested in hemi region
}; // MetaCovTest

class MetaDominantCovTest: public MetaCovTest{
 public:
  MetaDominantCovTest(int windowSize): MetaCovTest(windowSize) {
    this->modelName = "MetaDominantCov";
  }
  virtual int fit(DataConsolidator* dc) {
    dc->codeGenotypeForDominantModel(&geno);
    return fitWithGivenGenotype(geno, dc);
  }
 private:
  Matrix geno;
};

class MetaRecessiveCovTest: public MetaCovTest{
 public:
  MetaRecessiveCovTest(int windowSize): MetaCovTest(windowSize) {
    this->modelName = "MetaRecessiveCov";
  }
  virtual int fit(DataConsolidator* dc) {
    dc->codeGenotypeForRecessiveModel(&geno);
    return fitWithGivenGenotype(geno, dc);
  }
 private:
  Matrix geno;
};


#if 0
class MetaSkewTest: public ModelFitter{
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
  MetaSkewTest(int windowSize){
    this->modelName = "MetaSkew";
    this->numVariant = 0;
    this->nSample = -1;
    // this->mleVarY = -1.;
    this->fout = NULL;
    this->windowSize = windowSize;
    this->needToFitNullModel = true;
    result.addHeader("CHROM");
    result.addHeader("START_POS");
    result.addHeader("END_POS");
    result.addHeader("NUM_MARKER");
    result.addHeader("MARKER_POS");
    result.addHeader("SKEW");
  };
  ~MetaSkewTest(){
    while(queue.size() > 0 ) {
      if (isBinaryOutcome()) {
        printSkewForBinaryTrait(fout, queue);
      }
      queue.pop_front();
    }
  }

  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();
    Result& siteInfo = dc->getResult();

    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }
    if (genotype.rows == 0){
      fitOK = false;
      return -1;
    }
    this->useFamilyModel = dc->hasKinship();
    if (nSample < 0) { // unitialized
      // calculate variance of y
      nSample = genotype.rows;
      weight.Dimension(nSample);
    } else {
      if (nSample != genotype.rows){
        fprintf(stderr, "Sample size changed at [ %s:%s ]", siteInfo["CHROM"].c_str(), siteInfo["POS"].c_str());
        fitOK = false;
        return -1;
      }
    }

    // set weight
    if (useFamilyModel) {
      fitOK = false;
      return -1;
    }

    if (!isBinaryOutcome()) {
      static bool warningGiven = false;
      if (!warningGiven) {
        fprintf(stderr, "For quantative trait, it is not necessary to use MetaSkew model.\n");
        warningGiven = true;
      }
      return -1;
    } else { // binary case
      // fit null model
      if (this->needToFitNullModel || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) {
        copyCovariateAndIntercept(genotype.rows, covariate, &cov);
        copyPhenotype(phenotype, &this->pheno);
        fitOK = logistic.FitNullModel(cov, pheno, 100);
        if (!fitOK) return -1;
        needToFitNullModel = false;

        // skip store Z, as Z = this->cov
        // store V in weight
        for (int i = 0; i < nSample; ++i) {
          const double y = logistic.GetNullPredicted()[i];
          weight[i] = y * (1.0 - y) * (1.0 - 2.0 * y);
        }
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
      if (isBinaryOutcome()) {
        printSkewForBinaryTrait(fout, queue);
      }
      // if (isBinaryOutcome()) {
      //   // printCovarianceForBinaryTrait(fout, queue);
      // } else {
      //   printCovariance(fout, queue);
      // }
      queue.pop_front();
    };
    if (fitOK) {
      queue.push_back(loci);
      ++numVariant;
    }
    // result.writeValueLine(fp);
  };

 private:
  // check if @param g is polymorphic
  bool isVariant(const Genotype& g) {
    int n1 = g.size();
    if (n1 <= 1) return true;

    for(int i = 1; i < n1; ++i) {
      if (g[i] != g[0]) return true;
    }
    return false;
  }

  // get weighted 3rd moment
  double getMoment(const Genotype& g1, const Genotype& g2, const Genotype& g3) {
    double sum_ijk = 0.0;

    int n1 = g1.size();
    int n2 = g2.size();
    int n3 = g3.size();
    if (!(n1 == n2 && n2 == n3)){
      assert(false);
    }

    for(int i = 0; i < n1; ++i) {
      if (g1[i] == 0.0 || g2[i] == 0.0 || g3[i] == 0.0) continue;
      sum_ijk += g1[i] * g2[i] * g3[i] * weight[i];
    }
    return sum_ijk;
  }

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
   * print the skewness for the front of loci to the rest of loci
   */
  int printSkewForBinaryTrait(FileWriter* fp, const std::deque<Loci>& lociQueue){
    // skip monomorphic site
    if (!isVariant(lociQueue.front().geno)) {
      return 0;
    }

    // record polymorphic locations
    std::vector<std::deque<Loci>::const_iterator> polymorphicLoci;
    for (std::deque<Loci>::const_iterator iter = lociQueue.begin();
         iter != lociQueue.end();
         ++ iter) {
      if (isVariant(iter->geno)) {
        polymorphicLoci.push_back(iter);
      }
    }
    if (polymorphicLoci.empty()) return 0;

    // output results
    std::string s;
    std::vector<std::string> skew;

    bool hasSmallPvalue = true;
    Vector genoVec(lociQueue.front().geno.size());
    for(int i = 0; i < (int)lociQueue.front().geno.size(); ++i) {
      genoVec[i] = lociQueue.front().geno[i];
    }
    if (!logistic.TestCovariate(cov, pheno, genoVec)) {
      // fitting failed
      hasSmallPvalue = false;
    } else {
      hasSmallPvalue = (logistic.GetPvalue() < 0.1);
    }

    if (hasSmallPvalue) {
      // keep every combinations (currentPos, i, j)
      for (size_t i = 0; i != polymorphicLoci.size(); ++i) {
        for (size_t j = i; j != polymorphicLoci.size(); ++j) {
          const double val = getMoment(lociQueue.front().geno,
                                       (polymorphicLoci[i]->geno),
                                       (polymorphicLoci[j]->geno));
          if (val == 0.0) continue;

          s.clear();
          s += toString(i);
          s += ',';
          s += toString(j);
          s += ',';
          s += toString(val);

          skew.push_back(s);
        }
      }
    } else {
      // keep (currentPos, currentPos, i ) and
      // keep (currentPos, i, i) positions only
      for (size_t i = 0; i != polymorphicLoci.size(); ++i) {
        const double val = getMoment(lociQueue.front().geno,
                                     lociQueue.front().geno,
                                     (polymorphicLoci[i]->geno));
        if (val == 0.0) continue;

        s = "0,";
        s += toString(i);
        s += ',';
        s += toString(val);
        skew.push_back(s);
      }

      for (size_t i = 1; i != polymorphicLoci.size(); ++i) {
        const double val = getMoment(lociQueue.front().geno,
                                     (polymorphicLoci[i]->geno),
                                     (polymorphicLoci[i]->geno));
        if (val == 0.0) continue;

        s.clear();
        s += toString(i);
        s += ',';
        s += toString(i);
        s += ',';
        s += toString(val);

        skew.push_back(s);
      }
    }

    result.updateValue("CHROM", lociQueue.front().pos.chrom);
    result.updateValue("START_POS", lociQueue.front().pos.pos);
    result.updateValue("END_POS", lociQueue.back().pos.pos);
    /* fprintf(fp, "%d\t", idx); */
    result.updateValue("NUM_MARKER", (int) polymorphicLoci.size());

    s.clear();
    for(size_t i = 0; i != polymorphicLoci.size(); ++i) {
      if (i) s += ',';
      s+= toString(polymorphicLoci[i]->pos.pos);
    }
    result.updateValue("MARKER_POS", s);

    s.clear();
    stringJoin(skew, ':', &s);
    result.updateValue("SKEW", s);
    result.writeValueLine(fp);
    return 0;
  };

 private:
  std::deque< Loci> queue;
  int numVariant;
  int nSample;
  // Vector mleVarY;  // variance term of Y_i for i = 1,..., N th sample
  FileWriter* fout;
  int windowSize;
  Loci loci;
  bool fitOK;
  // Result result;
  bool useFamilyModel;
  Vector weight; // per individual weight
  LogisticRegressionScoreTest logistic;
  bool needToFitNullModel;
  Matrix cov; // covariate
  Vector pheno;
  std::map< std::pair<std::string, int>, bool> hasSmallPvalue;
}; // MetaSkewTest


class MetaKurtTest: public ModelFitter{
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
  MetaKurtTest(int windowSize){
    this->modelName = "MetaKurt";
    this->numVariant = 0;
    this->nSample = -1;
    // this->mleVarY = -1.;
    this->fout = NULL;
    this->windowSize = windowSize;
    this->mafThreshold = 0.05;
    result.addHeader("CHROM");
    result.addHeader("START_POS");
    result.addHeader("END_POS");
    result.addHeader("NUM_MARKER");
    result.addHeader("MARKER_POS");
    result.addHeader("KURT");
  };
  ~MetaKurtTest(){
    while(queue.size() > 0 ) {
      if (isBinaryOutcome()) {
        printKurtForBinaryTrait(fout, queue);
      }
      queue.pop_front();
    }
  }

  // fitting model
  int fit(DataConsolidator* dc) {
    // Matrix& phenotype = dc-> getPhenotype();
    Matrix& genotype = dc->getGenotype();
    // Matrix& covariate = dc->getCovariate();
    Result& siteInfo = dc->getResult();

    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }
    if (genotype.rows == 0){
      fitOK = false;
      return -1;
    }
    this->useFamilyModel = dc->hasKinship();
    if (nSample < 0) { // unitialized
      // calculate variance of y
      nSample = genotype.rows;
      weight.Dimension(nSample);
    } else {
      if (nSample != genotype.rows){
        fprintf(stderr, "Sample size changed at [ %s:%s ]", siteInfo["CHROM"].c_str(), siteInfo["POS"].c_str());
        fitOK = false;
        return -1;
      }
    }

    // set weight
    if (useFamilyModel) {
      fitOK = false;
      return -1;
    }

    if (!isBinaryOutcome()) {
      static bool warningGiven = false;
      if (!warningGiven) {
        fprintf(stderr, "For quantative trait, it is not necessary to use MetaKurt model.\n");
        warningGiven = true;
      }
      return -1;
    } else { // binary case
      /* // fit null model */
      /* if (this->needToFitNullModel || dc->isPhenotypeUpdated() || dc->isCovariateUpdated()) { */
      /*   copyCovariateAndIntercept(genotype.rows, covariate, &cov); */
      /*   copyPhenotype(phenotype, &this->pheno); */
      /*   fitOK = logistic.FitNullModel(cov, pheno, 100); */
      /*   if (!fitOK) return -1; */
      /*   needToFitNullModel = false; */

      /*   // skip store Z, as Z = this->cov */
      /*   // store V in weight */
      /*   for (int i = 0; i < nSample; ++i) { */
      /*     const double y = logistic.GetNullPredicted()[i]; */
      /*     weight[i] = y * (1.0 - y) * (1.0 - 2.0 * y); */
      /*   } */
      /* } */
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
      if (isBinaryOutcome()) {
        printKurtForBinaryTrait(fout, queue);
      }
      // if (isBinaryOutcome()) {
      //   // printCovarianceForBinaryTrait(fout, queue);
      // } else {
      //   printCovariance(fout, queue);
      // }
      queue.pop_front();
    };
    if (fitOK) {
      queue.push_back(loci);
      ++numVariant;
    }
    // result.writeValueLine(fp);
  };

 private:
  // check if MAF of @param g is larger than this->mafThreshold
  bool passMAFThreshold(const Genotype& g) {
    int n1 = g.size();
    if (n1 <= 0) return false;

    double s = 0.0;
    for(int i = 0; i < n1; ++i) {
      if (g[i]<0) continue;
      s += g[i];
    }
    double af = 0.5 * s / n1;
    if (af > .5) {
      af = 1.0 - af;
    }
    if (af < this->mafThreshold) return false;
    return true;
  }

  /** get weighted 4th moment
   * Take @param g1 and @param g2,
   * calculate \sum g1 * g1 * g2 * g2 to @param kurt1
   * calculate \sum g1 * g1 * g1 * g2 to @param kurt12
   */
  double getMoment(const Genotype& g1,
                   const Genotype& g2,
                   double* kurt1,
                   double* kurt2) {
    double sum_i2j2 = 0.0;
    double sum_i3j1 = 0.0;

    int n1 = g1.size();
    int n2 = g2.size();
    if (!(n1 == n2)){
      assert(false);
      return false;
    }

    for(int i = 0; i < n1; ++i) {
      if (g1[i] == 0.0 || g2[i] == 0.0) continue;
      double d = g1[i] * g1[i] * g2[i];
      sum_i2j2 += d * g1[i];
      sum_i3j1 += d * g2[i];
    }
    return true;
  }

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
   * print the kurtosis for the front of loci to the rest of loci
   */
  int printKurtForBinaryTrait(FileWriter* fp, const std::deque<Loci>& lociQueue){
    // skip monomorphic site
    if (!passMAFThreshold(lociQueue.front().geno)) {
      return 0;
    }

    // record polymorphic locations
    std::vector<std::deque<Loci>::const_iterator> polymorphicLoci;
    for (std::deque<Loci>::const_iterator iter = lociQueue.begin();
         iter != lociQueue.end();
         ++ iter) {
      if (this->passMAFThreshold(iter->geno)) {
        polymorphicLoci.push_back(iter);
      }
    }
    if (polymorphicLoci.empty()) return 0;

    // output results
    std::string s;
    std::vector<std::string> kurt;

    // keep every combinations (currentPos, i, j)
    double val_i2j2 = 0.0;
    double val_i3j1 = 0.0;
    for (size_t i = 0; i != polymorphicLoci.size(); ++i) {
      getMoment(lociQueue.front().geno,
                (polymorphicLoci[i]->geno),
                &val_i2j2,
                &val_i3j1);

      if (val_i2j2 == 0.0 && val_i3j1 == 0.0) continue;

      s.clear();
      s += toString(i);
      s += ',';
      s += toString(val_i2j2);
      s += ',';
      s += toString(val_i3j1);
      kurt.push_back(s);
    }

    result.updateValue("CHROM", lociQueue.front().pos.chrom);
    result.updateValue("START_POS", lociQueue.front().pos.pos);
    result.updateValue("END_POS", lociQueue.back().pos.pos);
    /* fprintf(fp, "%d\t", idx); */
    result.updateValue("NUM_MARKER", (int) polymorphicLoci.size());

    s.clear();
    for(size_t i = 0; i != polymorphicLoci.size(); ++i) {
      if (i) s += ',';
      s+= toString(polymorphicLoci[i]->pos.pos);
    }
    result.updateValue("MARKER_POS", s);

    s.clear();
    stringJoin(kurt, ':', &s);
    // fprintf(stderr, "s.size() = %zu\n", s.size());
    result.updateValue("KURT", s);
    result.writeValueLine(fp);
    return 0;
  };

 private:
  std::deque< Loci> queue;
  int numVariant;
  int nSample;
  // Vector mleVarY;  // variance term of Y_i for i = 1,..., N th sample
  FileWriter* fout;
  int windowSize;
  Loci loci;
  bool fitOK;
  // Result result;
  bool useFamilyModel;
  Vector weight; // per individual weight
  Matrix cov; // covariate
  Vector pheno;
  std::map< std::pair<std::string, int>, bool> hasSmallPvalue;
  double mafThreshold;
}; // MetaKurtTest

#endif

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
    // copy genotype column label
    for (int i = 0; i < genotype.cols; ++i){
      this->genotype.SetColumnLabel(i, genotype.GetColumnLabel(i));
    }
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
      if (strlen(genotype.GetColumnLabel(i)) == 0) {
        fprintf(fDump, "\tX%d", i);
      } else {
        fprintf(fDump, "\t%s", genotype.GetColumnLabel(i));
      }
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

#endif /* _MODELFITTER_H_ */
