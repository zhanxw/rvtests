#ifndef _MODEL_H_
#define _MODEL_H_

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <deque>

#include "libsrc/MathMatrix.h"

#include "base/Argument.h"
#include "base/ParRegion.h"
#include "base/RingMemoryPool.h"
#include "regression/MatrixRef.h"

#include "regression/BoltLMM.h"
#include "regression/EigenMatrixInterface.h"
#include "regression/FamSkat.h"
#include "regression/FastLMM.h"
#include "regression/FastMultipleTraitLinearRegressionScoreTest.h"
#include "regression/FirthRegression.h"
#include "regression/GrammarGamma.h"
#include "regression/LinearRegression.h"
#include "regression/LinearRegressionScoreTest.h"
#include "regression/LinearRegressionVT.h"
#include "regression/LogisticRegression.h"
#include "regression/LogisticRegressionScoreTest.h"
#include "regression/LogisticRegressionVT.h"
#include "regression/MatrixOperation.h"
#include "regression/MetaCov.h"
#include "regression/MultipleTraitLinearRegressionScoreTest.h"
#include "regression/MultivariateVT.h"
#include "regression/Skat.h"
#include "regression/SkatO.h"
#include "regression/Table2by2.h"
#include "regression/kbac_interface.h"

#include "src/DataConsolidator.h"
#include "src/LinearAlgebra.h"
#include "src/ModelFitter.h"
#include "src/ModelParser.h"
#include "src/ModelUtil.h"
#include "src/Permutation.h"
#include "src/Result.h"
#include "src/Summary.h"

#if 0
// may decrease speed.
#ifdef _OPENMP
#include <omp.h>
#pragma message "Enable multithread using OpenMP"
#endif
#endif
#include <omp.h>
DECLARE_BOOL_PARAMETER(hideCovar);

extern SummaryHeader* g_SummaryHeader;

typedef void (*CollapsingFunction)(Matrix& in, const std::vector<int>& idx,
                                   Matrix* out, int index);

// various collapsing method
// they all take people by marker matrix
// and they won't take special care of missing genotypes

// double getMarkerFrequency(Matrix& in, int col);
// void getMarkerFrequency(Matrix& in, std::vector<double>* freq);
double getMarkerFrequency(DataConsolidator* dc, int col);
void getMarkerFrequency(DataConsolidator* dc, std::vector<double>* freq);
double getMarkerFrequencyFromControl(Matrix& in, Vector& pheno, int col);

void cmcCollapse(DataConsolidator* dc, Matrix& in, Matrix* out);
void cmcCollapse(DataConsolidator* dc, Matrix& in, const std::vector<int>& idx,
                 Matrix* out, int index);

void zegginiCollapse(DataConsolidator* dc, Matrix& in, Matrix* out);
void zegginiCollapse(DataConsolidator* dc, Matrix& in,
                     const std::vector<int>& idx, Matrix* out, int index);

void fpCollapse(DataConsolidator* dc, Matrix& in, Matrix* out);

void madsonBrowningCollapse(DataConsolidator* dc, Matrix& genotype,
                            Vector& phenotype, Matrix* out);

void groupFrequency(const std::vector<double>& freq,
                    std::map<double, std::vector<int> >* group);
void convertToReferenceAlleleCount(Matrix& in, Matrix* g);

void makeVariableThreshodlGenotype(
    DataConsolidator* dc, Matrix& in, Matrix* out, std::vector<double>* freqOut,
    void (*collapseFunc)(DataConsolidator* dc, Matrix&, const std::vector<int>&,
                         Matrix*, int));

void appendHeritability(FileWriter* fp, const FastLMM& model);
void appendHeritability(FileWriter* fp, const GrammarGamma& model);

class SingleVariantWaldTest : public ModelFitter {
 public:
  SingleVariantWaldTest() : fitOK(false) {
    this->modelName = "SingleWald";
    result.addHeader("Test");
    result.addHeader("Beta");
    result.addHeader("SE");
    result.addHeader("Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& cov = dc->getCovariate();

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
    if (isMonomorphicMarker(genotype, 0)) {
      // Put monomorphic check after copying genotype to this->X
      // so that the output column "Test" can work
      fitOK = false;
      return -1;
    }

    if (!isBinaryOutcome()) {
      fitOK = linear.FitLinearModel(this->X, this->Y);
    } else {
      fitOK = logistic.FitLogisticModel(this->X, this->Y, 100);
    }
    return (fitOK ? 0 : 1);
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    // fprintf(fp, "Test\tBeta\tSE\tPvalue\n");
    result.writeHeaderLine(fp);
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    // skip interecept (column 0)
    for (int i = 1; i < this->X.cols; ++i) {
      if (FLAG_hideCovar && i > 1) {
        continue;
      }
      siteInfo.writeValueTab(fp);
      result.updateValue("Test", this->X.GetColumnLabel(i));
      if (fitOK) {
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
  }

 private:
  Matrix X;  // 1 + cov + geno
  Vector Y;  // phenotype
  LinearRegression linear;
  LogisticRegression logistic;
  bool fitOK;
};  // SingleVariantWaldTest

class SingleVariantFirthTest : public ModelFitter {
 public:
  SingleVariantFirthTest() : fitOK(false) {
    this->modelName = "SingleFirth";
    result.addHeader("Test");
    result.addHeader("Beta");
    result.addHeader("SE");
    result.addHeader("Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& cov = dc->getCovariate();

    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }
    if (!isBinaryOutcome()) {
      warnOnce(
          "Firth test does not support continuous outcomes. Results will be "
          "all NAs.");
      fitOK = false;
      return -1;
    }

    copyPhenotype(phenotype, &this->Y);
    if (cov.cols) {
      copyGenotypeWithCovariateAndIntercept(genotype, cov, &this->X);
    } else {
      copyGenotypeWithIntercept(genotype, &this->X);
    }
    if (isMonomorphicMarker(genotype, 0)) {
      // Put monomorphic check after copying genotype to this->X
      // so that the output column "Test" can work
      fitOK = false;
      return -1;
    }

    fitOK = firth.Fit(this->X, this->Y);
    // dumpToFile(X, "X");
    // dumpToFile(Y, "Y");
    return (fitOK ? 0 : 1);
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    // fprintf(fp, "Test\tBeta\tSE\tPvalue\n");
    result.writeHeaderLine(fp);
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    // skip interecept (column 0)
    for (int i = 1; i < this->X.cols; ++i) {
      siteInfo.writeValueTab(fp);
      result.clearValue();
      result.updateValue("Test", this->X.GetColumnLabel(i));
      if (fitOK) {
        result.updateValue("Beta", firth.GetCovEst()[i]);
        result.updateValue("SE", sqrt(firth.GetCovB()[i][i]));
        result.updateValue("Pvalue", firth.GetAsyPvalue()[i]);
      }
      result.writeValueLine(fp);
    }
  }

 private:
  Matrix X;  // 1 + cov + geno
  Vector Y;  // phenotype
  FirthRegression firth;
  bool fitOK;
};  // SingleVariantFirthTest

class SingleVariantScoreTest : public ModelFitter {
 public:
  SingleVariantScoreTest()
      : af(-1), nSample(-1), fitOK(false), needToFitNullModel(true) {
    this->modelName = "SingleScore";
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();

    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }
    nSample = genotype.rows;
    af = getMarkerFrequency(dc, 0);
    if (isMonomorphicMarker(genotype, 0)) {
      fitOK = false;
      return -1;
    }

    copyCovariateAndIntercept(genotype.rows, covariate, &cov);
    copyPhenotype(phenotype, &this->pheno);

    if (!isBinaryOutcome()) {
      if (needToFitNullModel || dc->isPhenotypeUpdated() ||
          dc->isCovariateUpdated()) {
        fitOK = linear.FitNullModel(cov, pheno);
        if (!fitOK) {
          warnOnce("Single variant score test failed in fitting null model.");
          return -1;
        }
        needToFitNullModel = false;
      }
      fitOK = linear.TestCovariate(cov, pheno, genotype);
    } else {
      if (needToFitNullModel || dc->isPhenotypeUpdated() ||
          dc->isCovariateUpdated()) {
        fitOK = logistic.FitNullModel(cov, pheno, 100);
        if (!fitOK) {
          warnOnce("Single variant score test failed in fitting null model.");
          return -1;
        }
        // calculateConstant(phenotype);
        needToFitNullModel = false;
      }
      fitOK = logistic.TestCovariate(cov, pheno, genotype);
    }
    return (fitOK ? 0 : -1);
  }

  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.addHeader("AF");
    result.addHeader("U");
    result.addHeader("V");
    result.addHeader("STAT");
    result.addHeader("DIRECTION");
    result.addHeader("EFFECT");
    result.addHeader("SE");
    result.addHeader("PVALUE");
    result.writeHeaderLine(fp);
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    result.clearValue();
    result.updateValue("AF", af);
    if (fitOK) {
      if (!isBinaryOutcome()) {
        const double u = linear.GetU()[0][0];
        const double v = linear.GetV()[0][0];
        result.updateValue("U", u);
        result.updateValue("V", v);
        result.updateValue("STAT", linear.GetStat());
        if (u != 0) {
          result.updateValue("DIRECTION", linear.GetU()[0][0] > 0 ? "+" : "-");
        }
        if (v > 0) {
          result.updateValue("EFFECT", linear.GetBeta()[0][0]);
          result.updateValue("SE", linear.GetSEBeta(0));
        }
        result.updateValue("PVALUE", linear.GetPvalue());
      } else {
        const double u = logistic.GetU()[0][0];
        const double v = logistic.GetV()[0][0];
        result.updateValue("U", u);
        result.updateValue("V", v);
        result.updateValue("STAT", logistic.GetStat());
        if (u != 0) {
          result.updateValue("DIRECTION",
                             logistic.GetU()[0][0] > 0 ? "+" : "-");
        }
        if (v > 0) {
          result.updateValue("EFFECT", u / v);
          result.updateValue("SE", 1.0 / sqrt(v));
        }
        result.updateValue("PVALUE", logistic.GetPvalue());
      }
    }
    result.writeValueLine(fp);
  }
  // don't need this
  // void calculateConstant(Matrix& phenotype);

 private:
  // double b;  // a constant
  double af;
  int nSample;
  Vector pheno;
  LinearRegressionScoreTest linear;
  LogisticRegressionScoreTest logistic;
  bool fitOK;
  bool needToFitNullModel;
  Matrix cov;
};  // SingleVariantScoreTest

class SingleVariantFisherExactTest : public ModelFitter {
 public:
  SingleVariantFisherExactTest() {
    this->modelName = "FisherExact";
    caseAC = caseAN = ctrlAC = ctrlAN = 0;
    fitOK = false;
  }
  virtual void countGenotype(int geno, int pheno) {
    if (geno == 0) {
      model.Increment(0, pheno);
      model.Increment(0, pheno);
    } else if (geno == 1) {
      model.Increment(0, pheno);
      model.Increment(1, pheno);
    } else if (geno == 2) {
      model.Increment(1, pheno);
      model.Increment(1, pheno);
    }
  }
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
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& cov = dc->getCovariate();

    if (genotype.cols != 1 || isMonomorphicMarker(genotype, 0)) {
      fitOK = false;
      return -1;
    }
    if (!isBinaryOutcome()) {
      warnOnce(
          "Fisher's exact test does not support continuous outcomes. Results "
          "will be all NAs.");
      fitOK = false;
      return -1;
    }
    if (cov.cols) {
      warnOnce(
          "Fisher's exact test does not support covariates. Results will be "
          "all NAs.");
      fitOK = false;
      return -1;
    }

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
      if (!(0 <= geno && geno <= 2)) continue;
      if (!(0 <= pheno && pheno <= 1)) continue;
      if (pheno == 1) {
        caseAC += geno;
        caseAN += 2;
      } else {
        ctrlAC += geno;
        ctrlAN += 2;
      }

      countGenotype(geno, pheno);
    }

    // step 2, calculate pvalue
    model.UpdateMarginSum();
    model.FullFastFisherExactTest();

    this->fitOK = true;
    return 0;
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (fitOK) {
      result.updateValue("N00", model.Get00());
      result.updateValue("N01", model.Get01());
      result.updateValue("N10", model.Get10());
      result.updateValue("N11", model.Get11());

      if (ctrlAN == 0) {
        result.updateValue("CtrlAF", 0);
      } else {
        result.updateValue("CtrlAF", 1.0 * ctrlAC / ctrlAN);
      }
      if (caseAN == 0) {
        result.updateValue("CaseAF", 0);
      } else {
        result.updateValue("CaseAF", 1.0 * caseAC / caseAN);
      }
      result.updateValue("PvalueTwoSide", model.getPExactTwoSided());
      result.updateValue("PvalueLess", model.getPExactOneSidedLess());
      result.updateValue("PvalueGreater", model.getPExactOneSidedGreater());
    }
    result.writeValueLine(fp);
  }
  void reset() {
    ModelFitter::reset();
    model.reset();
  }

 protected:
  Table2by2 model;

 private:
  int caseAC;
  int caseAN;
  int ctrlAC;
  int ctrlAN;

  bool fitOK;
};  // SingleVariantFisherExactTest

class SingleVariantDominantFisherExactTest
    : public SingleVariantFisherExactTest {
 public:
  SingleVariantDominantFisherExactTest() : SingleVariantFisherExactTest() {
    this->modelName = "DominantFisherExact";
  }

 public:
  virtual void countGenotype(int geno, int pheno) {
    if (geno == 0)
      model.Increment(0, pheno);
    else if (geno > 0)
      model.Increment(1, pheno);
  }
};

class SingleVariantFamilyScore : public ModelFitter {
 public:
  SingleVariantFamilyScore()
      : model(FastLMM::SCORE, FastLMM::MLE),
        fitOK(false),
        af(-1.),
        u(-1.),
        v(-1.),
        pvalue(-1.) {
    this->modelName = "FamScore";
    this->familyModel = true;
    result.addHeader("AF");
    result.addHeader("U.Stat");
    result.addHeader("V.Stat");
    result.addHeader("Pvalue");
    needToFitNullModel = true;
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    if (isBinaryOutcome()) {
      warnOnce(
          "Single variant score test (related individual) does not support "
          "binary outcomes. Results will be all NAs.");
      fitOK = false;
      return -1;
    }
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();

    if (genotype.cols != 1 || isMonomorphicMarker(genotype, 0)) {
      fitOK = false;
      return -1;
    }

    if (needToFitNullModel || dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK =
          (0 ==
                   model.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(),
                                      *dc->getKinshipSForAuto())
               ? true
               : false);
      if (!fitOK) return -1;
      needToFitNullModel = false;
    }

    fitOK = (0 ==
                     model.TestCovariate(cov, phenotype, genotype,
                                         *dc->getKinshipUForAuto(),
                                         *dc->getKinshipSForAuto())
                 ? true
                 : false);
    af = model.GetAF(*dc->getKinshipUForAuto(), *dc->getKinshipSForAuto(),
                     dc->getGenotype());
    u = model.GetUStat();
    v = model.GetVStat();
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
        result.updateValue("U.Stat", u);
        result.updateValue("V.Stat", v);
        result.updateValue("Pvalue", pvalue);
      }
    }
    result.writeValueLine(fp);
  }
  void writeFootnote(FileWriter* fp) {
    // appendHeritability(fp, model);
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
};  // end SingleVariantFamilyScore

class SingleVariantFamilyLRT : public ModelFitter {
 public:
  SingleVariantFamilyLRT()
      : model(FastLMM::LRT, FastLMM::MLE),
        fitOK(false),
        af(-1.),
        nullLogLik(-1.),
        altLogLik(-1.),
        pvalue(-1.) {
    this->modelName = "FamLRT";
    this->familyModel = true;
    result.addHeader("AF");
    result.addHeader("NullLogLik");
    result.addHeader("AltLogLik");
    result.addHeader("Pvalue");
    needToFitNullModel = true;
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    if (isBinaryOutcome()) {
      warnOnce(
          "Single variant likelihood ratio test (related individuals) does not "
          "support binary outcomes. Results will be all NAs.");
      fitOK = false;
      return -1;
    }
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();

    if (genotype.cols != 1 || isMonomorphicMarker(genotype, 0)) {
      fitOK = false;
      return -1;
    }

    if (needToFitNullModel || dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK =
          (0 ==
                   model.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(),
                                      *dc->getKinshipSForAuto())
               ? true
               : false);
      if (!fitOK) return -1;
      needToFitNullModel = false;
    }

    fitOK = (0 ==
             model.TestCovariate(cov, phenotype, genotype,
                                 *dc->getKinshipUForAuto(),
                                 *dc->getKinshipSForAuto()));
    af = model.GetAF(*dc->getKinshipUForAuto(), *dc->getKinshipSForAuto(),
                     dc->getGenotype());
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
    // appendHeritability(fp, model);
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
};  // end SingleVariantFamilyLRT

class SingleVariantFamilyGrammarGamma : public ModelFitter {
 public:
  SingleVariantFamilyGrammarGamma(GrammarGamma::AFMethod afMethod)
      : model(afMethod),
        needToFitNullModel(true),
        fitOK(false),
        af(-1.),
        beta(-1.),
        betaVar(-1.),
        pvalue(-1.) {
    this->modelName = "FamGrammarGamma";
    this->familyModel = true;
    result.addHeader("AF");
    result.addHeader("Beta");
    result.addHeader("BetaVar");
    result.addHeader("Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    if (isBinaryOutcome()) {
      warnOnce(
          "Single variant garmma-gamma test (related individuals) does not "
          "support binary outcomes. Results will be all NAs.");
      fitOK = false;
      return -1;
    }
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();

    if (genotype.cols != 1 || isMonomorphicMarker(genotype, 0)) {
      fitOK = false;
      return -1;
    }

    if (needToFitNullModel || dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK =
          (0 ==
                   model.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(),
                                      *dc->getKinshipSForAuto())
               ? true
               : false);
      if (!fitOK) return -1;
      needToFitNullModel = false;
    }

    fitOK = (0 ==
                     model.TestCovariate(cov, phenotype, genotype,
                                         *dc->getKinshipUForAuto(),
                                         *dc->getKinshipSForAuto())
                 ? true
                 : false);
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
    // appendHeritability(fp, model);
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
};  // SingleVariantFamilyGrammarGamma

class CMCTest : public ModelFitter {
 public:
  CMCTest() : fitOK(false), numVariant(-1) {
    this->modelName = "CMC";
    result.addHeader("NonRefSite");
#if 0
    result.addHeader("U.Stat");
    result.addHeader("V.Stat");
    result.addHeader("Effect");
    result.addHeader("SE");
#endif
    result.addHeader("Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
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
    for (int i = 0; i < phenotype.rows; i++) {
      pheno[i] = phenotype[i][0];
    }

    cmcCollapse(dc, genotype, &collapsedGenotype);

    if (isBinaryOutcome()) {
      fitOK = logistic.FitNullModel(cov, pheno, 100);
      if (!fitOK) {
        warnOnce("CMC test failed in fitting null model.");
        return -1;
      }
      fitOK = logistic.TestCovariate(cov, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    } else {
      fitOK = linear.FitNullModel(cov, pheno);
      if (!fitOK) {
        warnOnce("CMC test failed in fitting null model.");
        return -1;
      }
      fitOK = linear.TestCovariate(cov, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    }
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
  }

 private:
  /**
   * If the genotype is not exactly 0.0, we will count is as non-reference site
   */
  int totalNonRefSite() {
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
};  // CMCTest

class CMCWaldTest : public ModelFitter {
 public:
  CMCWaldTest() : fitOK(false), numVariant(-1) {
    this->modelName = "CMCWald";
    result.addHeader("NonRefSite");
    result.addHeader("Beta");
    result.addHeader("SE");
    result.addHeader("Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();

    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }

    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i < phenotype.rows; i++) {
      pheno[i] = phenotype[i][0];
    }

    cmcCollapse(dc, genotype, &collapsedGenotype);

    if (covariate.cols) {
      copyGenotypeWithCovariateAndIntercept(collapsedGenotype, covariate,
                                            &this->X);
    } else {
      copyGenotypeWithIntercept(collapsedGenotype, &this->X);
    }

    if (!isBinaryOutcome()) {
      fitOK = linear.FitLinearModel(X, phenotype);
    } else {
      fitOK = logistic.FitLogisticModel(X, phenotype, 100);
    }
    return (fitOK ? 0 : -1);
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.writeHeaderLine(fp);
  }
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
  }

 private:
  /**
   * If the genotype is not exactly 0.0, we will count is as non-reference site
   */
  int totalNonRefSite() {
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
};  // CMCWaldTest

class ZegginiWaldTest : public ModelFitter {
 public:
  ZegginiWaldTest() : fitOK(false), numVariant(-1) {
    this->modelName = "ZegginiWald";
    result.addHeader("Beta");
    result.addHeader("SE");
    result.addHeader("Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();

    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }

    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i < phenotype.rows; i++) {
      pheno[i] = phenotype[i][0];
    }

    zegginiCollapse(dc, genotype, &collapsedGenotype);

    if (covariate.cols) {
      copyGenotypeWithCovariateAndIntercept(collapsedGenotype, covariate,
                                            &this->X);
    } else {
      copyGenotypeWithIntercept(collapsedGenotype, &this->X);
    }

    if (!isBinaryOutcome()) {
      fitOK = linear.FitLinearModel(X, phenotype);
    } else {
      fitOK = logistic.FitLogisticModel(X, phenotype, 100);
    }
    return (fitOK ? 0 : -1);
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.writeHeaderLine(fp);
  }
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
  }

 private:
  Matrix collapsedGenotype;
  Matrix X;
  LogisticRegression logistic;
  LinearRegression linear;
  bool fitOK;
  int numVariant;
};  // ZegginiWaldTest

class CMCFisherExactTest : public ModelFitter {
 public:
  CMCFisherExactTest() : fitOK(false) {
    this->modelName = "CMCFisherExact";
    result.addHeader("N00");
    result.addHeader("N01");
    result.addHeader("N10");
    result.addHeader("N11");
    result.addHeader("PvalueTwoSide");
    result.addHeader("PvalueLess");
    result.addHeader("PvalueGreater");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& cov = dc->getCovariate();

    if (!isBinaryOutcome()) {
      warnOnce(
          "Fisher's exact test does not support continuous outcomes. Results "
          "will be all NAs.");
      fitOK = false;
      return -1;
    }
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    if (cov.cols) {
      warnOnce(
          "Fisher's exact test does not support covariates. Results will be "
          "all NAs.");
      fitOK = false;
      return -1;
    }

    // collapsing
    cmcCollapse(dc, genotype, &collapsedGenotype);

    // fit model
    // step 1, fit two by two table
    int numPeople = collapsedGenotype.rows;
    for (int i = 0; i < numPeople; i++) {
      int geno = collapsedGenotype[i][0];
      int pheno = phenotype[i][0];
      if (!(0 <= geno && geno <= 1)) continue;
      if (!(0 <= pheno && pheno <= 1)) continue;
      model.Increment(geno, pheno);
    }

    // step 2, calculate pvalue
    model.UpdateMarginSum();
    model.FullFastFisherExactTest();

    this->fitOK = true;
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
  }
  void reset() {
    ModelFitter::reset();
    model.reset();
  }

 private:
  Matrix collapsedGenotype;
  Table2by2 model;
  bool fitOK;
};  // CMCFisherExactTest

class ZegginiTest : public ModelFitter {
 public:
  ZegginiTest() : fitOK(false), numVariant(-1) {
    this->modelName = "Zeggini";
    result.addHeader("Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
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
    for (int i = 0; i < phenotype.rows; i++) {
      pheno[i] = phenotype[i][0];
    }

    zegginiCollapse(dc, genotype, &collapsedGenotype);

    if (isBinaryOutcome()) {
      fitOK = logistic.FitNullModel(cov, pheno, 100);
      if (!fitOK) {
        warnOnce("Zeggini test failed in fitting null model.");
        return -1;
      }
      fitOK = logistic.TestCovariate(cov, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    } else {
      fitOK = linear.FitNullModel(cov, pheno);
      if (!fitOK) {
        warnOnce("Zeggini test failed in fitting null model.");
        return -1;
      }
      fitOK = linear.TestCovariate(cov, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    }
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    // fprintf(fp, "Pvalue\n");
    result.writeHeaderLine(fp);
  }
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
  }

 private:
  Matrix collapsedGenotype;
  LogisticRegressionScoreTest logistic;
  LinearRegressionScoreTest linear;
  bool fitOK;
  int numVariant;
};  // ZegginiTest

class MadsonBrowningTest : public ModelFitter {
 public:
  MadsonBrowningTest(int nPerm, double alpha)
      : fitOK(false), numVariant(-1), perm(nPerm, alpha) {
    this->modelName = "MadsonBrowning";
    result.addHeader("Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();

    if (!isBinaryOutcome()) {
      warnOnce(
          "Madsen-Browning test does not support continuous outcomes. Results "
          "will be all NAs.");
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
    madsonBrowningCollapse(dc, genotype, pheno, &collapsedGenotype);

    fitOK = logistic.FitNullModel(cov, pheno, 100);
    if (!fitOK) {
      warnOnce("Madsen-Browning test failed in fitting null model.");
      return -1;
    }
    fitOK = logistic.TestCovariate(cov, pheno, collapsedGenotype);
    if (!fitOK) return -1;

    // record observed stat
    this->perm.init(logistic.GetStat());  // a chi-dist

    int failed = 0;
    while (this->perm.next()) {
      permute(&this->pheno);
      madsonBrowningCollapse(dc, genotype, pheno, &collapsedGenotype);
      fitOK = logistic.TestCovariate(collapsedGenotype, pheno);
      if (!fitOK) {
        if (failed < 10) {
          failed++;
          continue;
        } else {
          fitOK = false;
          return -1;
        }
      }
      // record new stats
      double pStat = logistic.GetStat();
      this->perm.add(pStat);
    }  // end permutation
    return (fitOK ? 0 : -1);
  }
  void reset() {
    ModelFitter::reset();
    this->perm.reset();
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);

    if (isBinaryOutcome()) {
      perm.writeHeaderTab(fp);
    } else {
      fp->write("Pvalue\n");
    }
    fp->write("\n");
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);

    if (isBinaryOutcome()) {
      perm.writeOutput(fp);
      fp->write("\n");
    } else {
      fp->write("NA\n");
    }
  }

 private:
  Matrix collapsedGenotype;
  Vector pheno;
  LogisticRegressionScoreTest logistic;
  bool fitOK;
  int numVariant;
  Permutation perm;
};  // MadsonBrowningTest

// Danyu Lin's method, using 1/sqrt(p(1-p)) as weight
// where p is estimated from all samples
class FpTest : public ModelFitter {
 public:
  FpTest() {
    this->modelName = "Fp";
    fitOK = false;
    numVariant = -1;
    result.addHeader("Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
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
    for (int i = 0; i < phenotype.rows; i++) {
      pheno[i] = phenotype[i][0];
    }

    fpCollapse(dc, genotype, &collapsedGenotype);

    if (isBinaryOutcome()) {
      fitOK = logistic.FitNullModel(cov, pheno, 100);
      if (!fitOK) {
        warnOnce("Fp test failed in fitting null model.");
        return -1;
      }
      fitOK = logistic.TestCovariate(cov, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    } else {
      fitOK = linear.FitNullModel(cov, pheno);
      if (!fitOK) {
        warnOnce("Fp test failed in fitting null model.");
        return -1;
      }
      fitOK = linear.TestCovariate(cov, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    }
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    // fprintf(fp, "Pvalue\n");
    result.writeHeaderLine(fp);
  }
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
  }

 private:
  Matrix collapsedGenotype;
  LogisticRegressionScoreTest logistic;
  LinearRegressionScoreTest linear;
  bool fitOK;
  int numVariant;
};  // FpTest

class RareCoverTest : public ModelFitter {
 public:
  RareCoverTest(int nPerm, double alpha)
      : fitOK(false), numVariant(-1), stat(-1), perm(nPerm, alpha) {
    this->modelName = "RareCover";
    this->result.addHeader("NumIncludeMarker");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();

    if (!isBinaryOutcome()) {
      warnOnce(
          "Rarecover test does not support continuous outcomes. Results will "
          "be all NAs.");
      fitOK = false;
      return -1;
    }
    if (covariate.cols != 0) {  // rare cover does not take covariate
      warnOnce(
          "Rarecover test does not support covariates. Results will be all "
          "NAs.");
      fitOK = false;
      return -1;
    }

    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    // use marker by people matrix for faster computation
    this->genotype.Transpose(genotype);
    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i < phenotype.rows; i++) {
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
    }
    fitOK = true;
    return (fitOK ? 0 : -1);
  }
  void reset() {
    ModelFitter::reset();
    this->perm.reset();
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    result.writeHeaderTab(fp);
    this->perm.writeHeaderLine(fp);
  }
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
  }
  /**
   * For a given genotype and phenotype, calculate RareCover stats, which
   * markers are selected
   * Here the genotype is: marker by people
   */
  double calculateStat(Matrix& genotype, Vector& phenotype,
                       std::set<int>* selectedIndex) {
    std::set<int>& selected = *selectedIndex;
    selected.clear();

    Vector c;  // collapsed genotype
    c.Dimension(phenotype.Length());
    c.Zero();

    double stat = -1;
    while ((int)selected.size() < genotype.rows) {
      int maxIdx = -1;
      double maxCorr = -1.0;
      double corr;
      for (int i = 0; i < genotype.rows; ++i) {
        if (selected.count(i)) continue;
        corr = calculateCorrelation(genotype[i], c, phenotype);
        if (corr > maxCorr) {
          maxCorr = corr;
          maxIdx = i;
        }
      }
      if (maxIdx < 0) {  // finish selection
        break;
      } else {
        if (maxCorr > stat) {
          // update selection
          stat = maxCorr;
          selected.insert(maxIdx);
          combine(&c, genotype[maxIdx]);
        } else {  // no select any new marker
          break;
        }
      }
    }
    return stat;
  }
  /**
   * Calculate correlatio of (g + collapsed, pheno)
   */
  double calculateCorrelation(Vector& g, Vector& collapsed, Vector& pheno) {
    double sum_g = 0.0;
    double sum_g2 = 0.0;
    double sum_p = 0.0;
    double sum_p2 = 0.0;
    double sum_gp = 0.0;
    int n = pheno.Length();
    for (int i = 0; i < n; ++i) {
      double geno = (g[i] + collapsed[i] > 0) ? 1.0 : 0.0;
      if (geno > 0.0) {
        sum_g += geno;
        sum_g2 += geno * geno;
        sum_gp += geno * pheno[i];
        sum_p += pheno[i];
        sum_p2 += pheno[i] * pheno[i];
      } else {
        sum_p += pheno[i];
        sum_p2 += pheno[i] * pheno[i];
      }
    }

    double cov_gp = sum_gp - sum_g * sum_p / n;
    double var_g = sum_g2 - sum_g * sum_g / n;
    double var_p = sum_p2 - sum_p * sum_p / n;
    double v = var_g * var_p;
    if (v < 1e-10) return 0.0;
    double corr = cov_gp / sqrt(v);
    return corr;
  }
  void combine(Vector* c, Vector& v) {
    int n = v.Length();
    for (int i = 0; i < n; ++i) {
      if ((*c)[i] + v[i] > 0) {
        (*c)[i] = 1.0;
      }
    }
  }

 private:
  Matrix genotype;
  std::set<int> selected;
  bool fitOK;
  int numVariant;
  double stat;
  Permutation perm;
};  // RareCoverTest

class CMATTest : public ModelFitter {
 public:
  CMATTest(int nPerm, double alpha) : perm(nPerm, alpha) {
    this->modelName = "CMAT";
    N_A = N_U = m_A = m_U = M_A = M_U = 0;
    fitOK = false;
    stat = -1.;
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& covariate = dc->getCovariate();

    if (!isBinaryOutcome()) {
      warnOnce(
          "CMAT test does not support continuous outcomes. Results will be all "
          "NAs.");
      fitOK = false;
      return -1;
    }
    if (covariate.cols != 0) {
      warnOnce(
          "CMAT exact test does not support covariates. Results will be all "
          "NAs.");
      fitOK = false;
      return -1;
    }

    // we use equal weight
    this->stat = this->calculateStat(dc, &N_A, &N_U, &m_A, &m_U, &M_A, &M_U);
    this->perm.init(this->stat);

    // permutation part
    double d1, d2, d3, d4, d5, d6;  // just used in permutation
    while (this->perm.next()) {
      permute(&pheno);
      // record new stats
      double pStat = this->calculateStat(dc, &d1, &d2, &d3, &d4, &d5, &d6);
      this->perm.add(pStat);
    }  // end permutation
    fitOK = true;
    return (fitOK ? 0 : -1);
  }
  void reset() {
    ModelFitter::reset();
    this->perm.reset();
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    if (isBinaryOutcome()) {  /// cmat only takes binary output
      this->perm.writeHeaderLine(fp);
    }
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (isBinaryOutcome()) {
      this->perm.writeOutputLine(fp);
    }
  }
  double calculateStat(DataConsolidator* dc, double* p_N_A, double* p_N_U,
                       double* p_m_A, double* p_m_U, double* p_M_A,
                       double* p_M_U) {
    Matrix& phenotype = dc->getPhenotype();
    copyPhenotype(phenotype, &this->pheno);
    Matrix& genotype = dc->getGenotype();

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

    for (int i = 0; i < pheno.Length(); ++i) {
      if (phenotype[i] == 1) {
        ++N_A;
      } else {
        ++N_U;
      }
    }
    for (int i = 0; i < genotype.cols; ++i) {
      // for each marker, get its allele frequency
      double af = getMarkerFrequency(dc, i);
      bool flip = false;
      if (af > 0.5) {
        flip = true;
      }
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
      }
    }
    if (N_A + N_U == 0.0) return 0.0;
    int numMarker = genotype.cols;
    return (N_A + N_U) / (2 * N_A * N_U * numMarker) * (m_A * M_U - m_U * M_A) *
           (m_A * M_U - m_U * M_A) / (m_A + m_U) / (M_A + M_U);
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
};  // CMATTest

#if 0
class UStatTest{
};// UStatTest
#endif

/**
 * Implementation of Alkes Price's VT with modifications
 *
 * 1. The original paper use reference allele counts and allele frequency
 *threshold.
 *    We use minor allele frequency and minor allele frequency threshold.
 * 2. The original permutation test p value = (x + 1) / (P + 1) and x is the
 *frequency
 *    of larger z from permutation, so it's a one-sided test.
 *    We use absolution value |z|, and we use a two-sided test.
 * 3. Original paper use 1000 permutation, but here we increase this number.
 * 4. Original paper and here both treat binary pheontype as continuous variable
 */
class VariableThresholdPrice : public ModelFitter {
 public:
  VariableThresholdPrice(int nPerm, double alpha)
      : fitOK(false), zmax(-1.), optimalFreq(-1), perm(nPerm, alpha) {
    this->modelName = "VariableThresholdPrice";
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();
    Vector& weight = dc->getWeight();

    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    if (covariate.cols != 0) {
      warnOnce(
          "Price's Variable Threshold method does not take covariate. Results "
          "will be all NAs");
    }
    if (weight.Length() != 0) {
      warnOnce(
          "Price's Variable Threshold method cannot take weights at this "
          "moment. Unweighted results will be calculated.");
    }

    // calculate allele frequency
    makeVariableThreshodlGenotype(dc, genotype, &this->sortedBurden,
                                  &this->freq, zegginiCollapse);
    transposeInPlace(
        &this->sortedBurden);  // now each row is a collapsed genoype at
    // certain frequency cutoff
    copyPhenotype(phenotype, &this->phenotype);

    this->zmax = -999.0;
    centerVector(&this->phenotype);
    // use unpermutated data
    if (calculateZ(this->phenotype, this->sortedBurden, weight, &this->zmax,
                   &this->optimalFreq)) {
      fitOK = false;
      return -1;
    }

    // begin permutation
    this->perm.init(fabs(zmax));
    double zPerm = 9999;
    double optFreqPerm;

    while (this->perm.next()) {
      permute(&this->phenotype);
      if (calculateZ(this->phenotype, this->sortedBurden, weight, &zPerm,
                     &optFreqPerm))
        continue;
      this->perm.add(fabs(zPerm));
    }

    fitOK = true;
    return 0;
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    fp->write("\tOptFreq\tZmax\t");
    this->perm.writeHeader(fp);
    fp->write("\n");
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    fp->printf("\t%g\t%g\t", this->optimalFreq, this->zmax);
    this->perm.writeOutputLine(fp);
    // fprintf(fp, "\n");
  }
  void reset() {
    fitOK = false;
    ModelFitter::reset();
    this->perm.reset();
  }

 private:
  void transposeInPlace(Matrix* g) {
    Matrix tmp = *g;
    g->Transpose(tmp);
  }

  int calculateZ(Vector& pheno, Matrix& sortedBurden, Vector& weight,
                 double* zmax, double* optimFreq) {
    assert(pheno.Length() == sortedBurden.cols);
    assert(sortedBurden.rows > 1);
    assert(weight.Length() == 0);
    assert(zmax);
    assert(optimFreq);
    assert((int)this->freq.size() == sortedBurden.rows);

    for (int i = 0; i < sortedBurden.rows; ++i) {
      const double z =
          calculateZthreshold(pheno, this->sortedBurden[i], weight);
      if (fabs(z) > *zmax || i == 0) {
        *zmax = fabs(z);
        *optimFreq = this->freq[i];
      }
    }
    return 0;
  }

  double calculateZthreshold(Vector& y, Vector& x, Vector& weight) {
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
    double sd = sqrt(getVariance(x));
    if (sd != 0) {
      ret /= sd;
    }
    return ret;
  }

 private:
  Matrix sortedBurden;
  std::vector<double> freq;

  bool fitOK;
  Vector phenotype;
  double zmax;
  double optimalFreq;  // the frequency cutoff in unpermutated data which give
  // smallest pvalue
  Permutation perm;
};  // VariableThresholdPrice

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
  }
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
  }
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
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    // this->result = siteInfo;
    // result.writeHeaderTab(fp);
    siteInfo.writeHeaderTab(fp);
    model[0].writeHeader(fp, result);
  }
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
  }
  void reset() {
    fitOK = true;
    for (int i = 0; i < this->modelLen; i++)
      model[i].reset();
  }
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

class VTCMC : public ModelFitter {
 public:
  VTCMC() : fitOK(false), needToFitNullModel(true) {
    this->modelName = "VTCMC";
    result.addHeader("MAF");
    result.addHeader("U");
    result.addHeader("V");
    result.addHeader("OptimMAF");
    result.addHeader("Effect");
    result.addHeader("Pvalue");
  }
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();

    copy(phenotype, &this->pheno);

    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    // copyPhenotype(phenotype, &this->phenotype);
    int ret = fitNullModel(dc);
    if (ret) {
      warnOnce("VT CMC test failed in fitting null model.");
      return -1;
    }
    // convertToMinorAlleleCount(genotype, &geno);
    // freq.clear();
    // for (int i = 0; i < genotype.rows; ++i){
    //   freq.push_back(getMarkerFrequency(geno, i));
    // }
    // groupFrequency(freq, &freqGroup);
    freq.clear();
    makeVariableThreshodlGenotype(dc, geno, &sortedGenotype, &freq,
                                  cmcCollapse);

    // rearrangeGenotypeByFrequency(genotype, &sortedGenotype, &this->freq);
    if (!isBinaryOutcome()) {
      ret = linear.TestCovariate(covariate, pheno, sortedGenotype);
    } else {
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
  }
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
    result.updateValue("MAF", floatToString(freq));
    result.updateValue("U", floatToString(U));
    result.updateValue("V", floatToString(V));
    result.updateValue("OptimMAF", optimFreq);
    result.updateValue("Effect", effect);
    result.updateValue("Pvalue", pValue);

    result.writeValueLine(fp);
  }
  void reset() {
    ModelFitter::reset();
    fitOK = false;
  }

 private:
  int fitNullModel(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();
    if (this->needToFitNullModel || dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
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
class AnalyticVT : public ModelFitter {
 public:
  typedef enum { UNRELATED = 0, RELATED = 1 } Type;

 public:
  AnalyticVT(AnalyticVT::Type type)
      : lmm(FastLMM::SCORE, FastLMM::MLE),
        fitOK(false),
        needToFitNullModel(true) {
    this->type = type;
    if (type == UNRELATED) {
      this->modelName = "AnalyticVT";
      this->familyModel = false;
    } else {
      this->modelName = "FamAnalyticVT";
      this->familyModel = true;
    }

    result.addHeader("MinMAF");
    result.addHeader("MaxMAF");
    result.addHeader("OptimMAF");
    result.addHeader("OptimNumVar");
    result.addHeader("U");
    result.addHeader("V");
    result.addHeader("Stat");
    result.addHeader("Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);

    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }

    if (isBinaryOutcome()) {
      // does not support binary outcomes
      warnOnce(
          "Analytic VT test does not support binary outcomes. Results will be "
          "all NAs.");
      return -1;
    }

    // obtain frequency, u and v per variant
    // const int nSample = genotype.rows;
    const int nVariant = genotype.cols;
    this->af.Dimension(nVariant);

    this->useFamilyModel = dc->hasKinship();
    if (this->useFamilyModel ^ (this->type == RELATED)) {
      // model and data does not match
      warnOnce("Analytic VT test has internal error!");
      return -1;
    }

    if (!this->useFamilyModel) {
      // calculate af
      for (int i = 0; i < nVariant; ++i) {
        af[i] = getMarkerFrequency(dc, i);
      }

      // adjust covariates
      if (needToFitNullModel || dc->isPhenotypeUpdated() ||
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
      if (needToFitNullModel || dc->isPhenotypeUpdated() ||
          dc->isCovariateUpdated()) {
        fitOK = lmm.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(),
                                 *dc->getKinshipSForAuto()) == 0;
        if (!fitOK) {
          warnOnce("Analytic VT test failed in fitting null model (LMM).");
          return -1;
        }
        needToFitNullModel = false;
      }

      for (int i = 0; i < nVariant; ++i) {
        af[i] = lmm.FastGetAF(*dc->getKinshipUForAuto(),
                              *dc->getKinshipSForAuto(), genotype, i);
        // fprintf(stderr, "af[%d] = %g\n", i, af[i]);
      }
      lmm.CalculateUandV(cov, phenotype, genotype, *dc->getKinshipUForAuto(),
                         *dc->getKinshipSForAuto(), &u, &v);
    }
    if (mvvt.compute(af, u, v)) {
      warnOnce("Analytic VT test failed in computing multivariate statistics.");
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
};  // AnalyticVT

class FamCMC : public ModelFitter {
 public:
  FamCMC()
      : needToFitNullModel(true),
        numVariant(0),
        u(-1.),
        v(-1.),
        af(-1.),
        effect(-1.),
        pvalue(-1.),
        lmm(FastLMM::SCORE, FastLMM::MLE),
        fitOK(false) {
    this->modelName = "FamCMC";
    this->familyModel = true;
    result.addHeader("NumSite");
    result.addHeader("AF");
    result.addHeader("U");
    result.addHeader("V");
    result.addHeader("Effect");
    result.addHeader("Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    if (isBinaryOutcome()) {
      warnOnce(
          "CMC test (for related individuals) does not support binary "
          "outcomes. Results will be "
          "all NAs.");
      fitOK = false;
      return -1;
    }

    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();

    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }

    if (needToFitNullModel || dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK =
          (0 ==
                   lmm.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(),
                                    *dc->getKinshipSForAuto())
               ? true
               : false);
      if (!fitOK) {
        warnOnce(
            "CMC test (for related individuals) failed in fitting null model "
            "(LMM).");
        return -1;
      }
      needToFitNullModel = false;
    }

    cmcCollapse(dc, genotype, &collapsedGenotype);

    // dumpToFile(genotype, "genotype");
    // dumpToFile(collapsedGenotype, "collapsedGenotype");

    fitOK = (0 ==
             lmm.TestCovariate(cov, phenotype, collapsedGenotype,
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
    // appendHeritability(fp, lmm);
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

class FamZeggini : public ModelFitter {
 public:
  FamZeggini()
      : needToFitNullModel(true),
        numVariant(-1),
        u(-1.),
        v(-1.),
        af(-1.),
        effect(-1.),
        pvalue(-1.),
        lmm(FastLMM::SCORE, FastLMM::MLE),
        fitOK(false) {
    this->modelName = "FamZeggini";
    this->familyModel = true;
    result.addHeader("NumSite");
    result.addHeader("MeanBurden");
    result.addHeader("U");
    result.addHeader("V");
    result.addHeader("Effect");
    result.addHeader("Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    if (isBinaryOutcome()) {
      warnOnce(
          "Zeggini test (for related individuals) does not support binary "
          "outcomes. Results will be "
          "all NAs.");
      fitOK = false;
      return -1;
    }

    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();

    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }

    if (needToFitNullModel || dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK =
          (0 ==
                   lmm.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(),
                                    *dc->getKinshipSForAuto())
               ? true
               : false);
      if (!fitOK) {
        warnOnce(
            "Zeggini test (for related individuals) failed in fitting null "
            "model (LMM).");
        return -1;
      }
      needToFitNullModel = false;
    }

    zegginiCollapse(dc, genotype, &collapsedGenotype);

    fitOK = (0 ==
             lmm.TestCovariate(cov, phenotype, collapsedGenotype,
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
    // appendHeritability(fp, lmm);
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

class FamFp : public ModelFitter {
 public:
  FamFp()
      : needToFitNullModel(true),
        numVariant(0),
        u(-1.),
        v(-1.),
        af(-1.),
        effect(-1.),
        pvalue(-1.),
        lmm(FastLMM::SCORE, FastLMM::MLE),
        fitOK(false) {
    this->modelName = "FamFp";
    this->familyModel = true;
    result.addHeader("NumSite");
    result.addHeader("AF");
    result.addHeader("U");
    result.addHeader("V");
    result.addHeader("Effect");
    result.addHeader("Pvalue");
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    if (isBinaryOutcome()) {
      warnOnce(
          "Fp test (for related individuals) does not support binary "
          "outcomes. Results will be "
          "all NAs.");
      fitOK = false;
      return -1;
    }

    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();

    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }

    if (needToFitNullModel || dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK =
          (0 ==
                   lmm.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(),
                                    *dc->getKinshipSForAuto())
               ? true
               : false);
      if (!fitOK) {
        warnOnce(
            "Fp test (for related individuals) failed in fitting null model "
            "(LMM).");
        return -1;
      }
      needToFitNullModel = false;
    }

    fpCollapse(dc, genotype, &collapsedGenotype);

    // dumpToFile(genotype, "genotype");
    // dumpToFile(collapsedGenotype, "collapsedGenotype");

    fitOK = (0 ==
             lmm.TestCovariate(cov, phenotype, collapsedGenotype,
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
    // appendHeritability(fp, lmm);
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

class SkatTest : public ModelFitter {
 public:
  /* SkatTest(const std::vector<std::string>& param) { */
  SkatTest(int nPerm, double alpha, double beta1, double beta2)
      : fitOK(false), pValue(-1.), stat(-1.), perm(nPerm, alpha) {
    this->usePermutation = nPerm > 0;
    this->beta1 = beta1;
    this->beta2 = beta2;
    this->modelName = "Skat";
    this->needToFitNullModel = true;
  }
  void reset() {
    ModelFitter::reset();
    this->skat.Reset();
    this->perm.reset();
    stat = -9999;
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();
    // not use dc->getWeight(), but use model specific weight
    // Vector& weight = dc->getWeight();

    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    // fill it weight
    // NOTE: our frequency calculation is slightly different than SKAT, so we
    // will need to adjust it back
    weight.Dimension(genotype.cols);
    for (int i = 0; i < weight.Length(); i++) {
      double freq = getMarkerFrequency(dc, i);
      // fprintf(stderr, "freq[%d] = %g\n", i, freq);
      if (freq > 0.5) {  // convert to MAF
        freq = 1.0 - freq;
      }
      if (freq > 1e-30) {  // avoid dividing zero
        weight[i] = gsl_ran_beta_pdf(
            freq, this->beta1,
            this->beta2);  /// default SKAT use beta(MAF, 1, 25)
        weight[i] *= weight[i];
        // fprintf(stderr, "weight(%d, %d, %f ) = %f\n", 1, 25, freq,
        // weight[i]);
      } else {
        weight[i] = 0.0;
      }
    }

    Vector phenoVec;
    copyPhenotype(phenotype, &phenoVec);

    // ynull is mean of y (removing genotypes) in model Ynull ~ X (aka Ynull ~ X
    // + 0.0 * G )
    Matrix cov;
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);

    // this part caluclate null model, only need to do once
    if (needToFitNullModel || dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
      if (isBinaryOutcome()) {
        fitOK = logistic.FitLogisticModel(cov, phenoVec, 100);
        if (!fitOK) {
          warnOnce("SKAT test failed in fitting null model (logistic model).");
          return -1;
        }
        ynull = logistic.GetPredicted();
        v = logistic.GetVariance();
      } else {
        fitOK = linear.FitLinearModel(cov, phenoVec);
        if (!fitOK) {
          warnOnce("SKAT test failed in fitting null model (linear model).");
          return -1;
        }
        ynull = linear.GetPredicted();
        v.Dimension(genotype.rows);
        for (int i = 0; i < genotype.rows; ++i) {
          v[i] = linear.GetSigma2();
        }
      }
      this->res.Dimension(ynull.Length());
      for (int i = 0; i < this->ynull.Length(); ++i) {
        this->res[i] = phenoVec[i] - this->ynull[i];
      }
      needToFitNullModel = false;
    }

    // get Pvalue
    skat.Fit(res, v, cov, genotype, weight);
    this->stat = skat.GetQ();
    this->pValue = skat.GetPvalue();

    // permuation part
    if (this->usePermutation) {
      permutedRes = res;
      this->perm.init(this->stat);

      double s;
      while (this->perm.next()) {
        permute(&permutedRes);
        s = skat.GetQFromNewResidual(permutedRes);
        this->perm.add(s);
      }
    }
    fitOK = true;
    return 0;
  }
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
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (!fitOK) {
      fp->write("NA\tNA");
      if (usePermutation) {
        fp->write("\tNA\tNA\tNA\tNA\tNA\tNA");
      }
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
  }

 private:
  bool needToFitNullModel;
  double beta1;
  double beta2;
  // Matrix X; // n by (p+1) matrix, people by covariate (note intercept is
  // needed);
  Vector v;
  Vector weight;
  LogisticRegression logistic;
  LinearRegression linear;
  Vector ynull;
  Vector res;          // residual under the null
  Vector permutedRes;  // residual under the null
  Skat skat;
  bool fitOK;
  double pValue;

  bool usePermutation;
  double stat;
  Permutation perm;
};  // SkatTest

class SkatOTest : public ModelFitter {
 public:
  SkatOTest(double beta1, double beta2) : fitOK(false) {
    this->beta1 = beta1;
    this->beta2 = beta2;
    this->modelName = "SkatO";
    this->needToFitNullModel = true;
  }
  void reset() {
    ModelFitter::reset();
    this->skato.Reset();
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();
    // not use dc->getWeight(), but use model specific weight
    // Vector& weight = dc->getWeight();

    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    // fill it weight
    weight.Dimension(genotype.cols);
    for (int i = 0; i < weight.Length(); i++) {
      double freq = getMarkerFrequency(dc, i);
      // fprintf(stderr, "freq[%d] = %g\n", i, freq);
      if (freq > 0.5) {  // convert to MAF
        freq = 1.0 - freq;
      }
      if (freq > 1e-30) {  // avoid dividing zero
        weight[i] = gsl_ran_beta_pdf(
            freq, this->beta1,
            this->beta2);  /// default SKAT use beta(MAF, 1, 25)
      } else {
        weight[i] = 0.0;
      }
    }

    Vector phenoVec;
    copyPhenotype(phenotype, &phenoVec);

    // ynull is the mean of y (removing genotypes) in the model:
    // Ynull ~ X (aka Ynull ~ X + 0.0 * G )
    Matrix cov;
    copyCovariateAndIntercept(genotype.rows, covariate, &cov);

    // this part caluclate null model, only need to do once
    if (needToFitNullModel || dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
      if (isBinaryOutcome()) {
        fitOK = logistic.FitLogisticModel(cov, phenoVec, 100);
        if (!fitOK) {
          warnOnce("SKAT test failed in fitting null model (logistic model).");
          return -1;
        }
        ynull = logistic.GetPredicted();
        v = logistic.GetVariance();
      } else {
        fitOK = linear.FitLinearModel(cov, phenoVec);
        if (!fitOK) {
          warnOnce("SKAT test failed in fitting null model (linear model).");
          return -1;
        }
        ynull = linear.GetPredicted();
        v.Dimension(genotype.rows);
        for (int i = 0; i < genotype.rows; ++i) {
          v[i] = linear.GetSigma2();
        }
      }
      this->res.Dimension(ynull.Length());
      for (int i = 0; i < this->ynull.Length(); ++i) {
        this->res[i] = phenoVec[i] - this->ynull[i];
      }
      needToFitNullModel = false;
    }

    // perform calculation
    if (!isBinaryOutcome()) {
      fitOK = skato.Fit(res, v, cov, genotype, weight, "C") == 0;
    } else {
      fitOK = skato.Fit(res, v, cov, genotype, weight, "D") == 0;
    }
    return 0;
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    fp->write("Q\trho\tPvalue\n");
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (!fitOK) {
      fp->write("NA\tNA\tNA\n");
    } else {
      fp->printf("%g\t%g\t%g\n", this->skato.GetQ(), this->skato.GetRho(),
                 this->skato.GetPvalue());
    }
  }

 private:
  bool needToFitNullModel;
  double beta1;
  double beta2;
  Vector v;
  Vector weight;
  LogisticRegression logistic;
  LinearRegression linear;
  Vector ynull;
  Vector res;  // residual under the null
  SkatO skato;
  bool fitOK;
};  // SkatOTest

class KBACTest : public ModelFitter {
 public:
  KBACTest(int nPerm, double alpha)
      : nPerm(nPerm),
        alpha(alpha),
        xdatIn(NULL),
        ydatIn(NULL),
        mafIn(NULL),
        xcol(0),
        ylen(0),
        nn(0),
        qq(0),
        aa(alpha),
        mafUpper(1.),
        twosided(1),
        fitOK(false),
        pValue(-1.) {
    this->modelName = "Kbac";
  }
  ~KBACTest() {
    if (this->xdatIn) delete[] this->xdatIn;
    if (this->ydatIn) delete[] this->ydatIn;
    if (this->mafIn) delete[] this->mafIn;
  }
  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    fp->write("Pvalue\n");
  }
  void reset() {
    // clear_kbac_test();
    ModelFitter::reset();
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();

    if (!isBinaryOutcome()) {
      warnOnce(
          "KBAC test does not support continuous outcomes. Results will be "
          "all NAs.");
      fitOK = false;
      return -1;
    }
    if (covariate.cols != 0) {
      warnOnce(
          "KBAC test does not support covariates. Results will be all NAs.");
      fitOK = false;
      return -1;
    }

    this->resize(genotype.rows, genotype.cols);
    this->nn = this->nPerm;
    this->qq = 1;
    this->aa = this->alpha;
    this->mafUpper = 1.0;  // no need to further prune alleles
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
        /*   fprintf(stderr, "i=%d, j=%d, genotype=%g\n", i,j,genotype[i][j]);
         */
        /*   fprintf(stderr, "j*genotype.rows +i = %d, xdatIn[..] = %g\n", j *
         * genotype.rows + i, xdatIn[j * genotype.rows + i]); */
        /* } */
      }
    }
    for (int i = 0; i < genotype.rows; ++i) {
      ydatIn[i] = phenotype[i][0];
    }
    for (int j = 0; j < genotype.cols; ++j) {
      mafIn[j] = getMarkerFrequency(dc, j);
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
    set_up_kbac_test(&nn, &qq, &aa, &mafUpper, xdatIn, ydatIn, mafIn, &xcol,
                     &ylen);
    /**
     * pvalue: results are here
     * twosided: two sided tests or one sided test
     */
    this->pValue =
        9.0;  // it is required to pass pvalue = 9.0 (see interface.R)
    do_kbac_test(&this->pValue, &this->twosided);
    // fprintf(stderr, "pValue = %g, twoside = %d\n", this->pValue,
    // this->twosided);
    clear_kbac_test();
    this->fitOK = true;
    return 0;
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeValueTab(fp);
    if (!fitOK) {
      fp->write("NA\n");
    } else {
      fp->printf("%f\n", this->pValue);
    }
  }
  void resize(int numPeople, int numMarker) {
    bool resized = false;
    if (numPeople != this->ylen) {
      delete[] this->ydatIn;
      this->ydatIn = new double[numPeople];
      this->ylen = numPeople;
      resized = true;
    }
    if (numMarker != this->xcol) {
      delete[] this->mafIn;
      this->mafIn = new double[numMarker];
      this->xcol = numMarker;
      resized = true;
    }

    if (resized) {
      delete[] this->xdatIn;
      this->xdatIn = new double[numPeople * numMarker];
    }
  }

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
};  // KBACTest

class FamSkatTest : public ModelFitter {
 public:
  FamSkatTest(double beta1, double beta2)
      : needToFitNullModel(true), fitOK(false), pValue(-1.), stat(-1.) {
    this->beta1 = beta1;
    this->beta2 = beta2;
    this->modelName = "FamSkat";
    this->familyModel = true;
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getFlippedToMinorPolymorphicGenotype();
    Matrix& covariate = dc->getCovariate();
    // not use dc->getWeight(), but use model specific weight
    // Vector& weight = dc->getWeight();

    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }

    if (isBinaryOutcome()) {
      warnOnce(
          "SKAT test (for related individuals) does not support binary "
          "outcomes. Results will be "
          "all NAs.");
      fitOK = false;
      return -1;
    }
    const bool useFamilyModel = dc->hasKinship();
    if (!useFamilyModel) {
      warnOnce(
          "SKAT test (for related individuals) cannot find kinship. Results "
          "will be "
          "all NAs.");
      fitOK = false;
      return -1;
    }

    // adjust covariates
    if (needToFitNullModel || dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK = skat.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(),
                                *dc->getKinshipSForAuto()) == 0;
      if (!fitOK) {
        warnOnce(
            "SKAT test (for related individuals) failed in fitting null model "
            "(SKAT)");
        return -1;
      }
      needToFitNullModel = false;
    }

    // get Pvalue
    fitOK = skat.TestCovariate(cov, phenotype, genotype, weight,
                               *dc->getKinshipUForAuto(),
                               *dc->getKinshipSForAuto()) == 0;
    if (!fitOK) {
      return -1;
    }

    this->stat = skat.GetQ();
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
    if (!fitOK) {
      fp->write("NA\tNA");
      fp->write("\n");
    } else {
      // binary outcome and quantative trait are similar output
      fp->printf("%g\t%g", this->skat.GetQ(), this->pValue);
      fp->write("\n");
    }
  }

 private:
  Vector weight;
  bool needToFitNullModel;
  double beta1;
  double beta2;
  FamSkat skat;
  bool fitOK;
  double pValue;
  Matrix cov;
  double stat;
};  // FamSkatTest

//////////////////////////////////////////////////////////////////////
// Meta-analysis based methods
//////////////////////////////////////////////////////////////////////

// output files for meta-analysis
// this class handle 8 cases:
// binary/qtl * unrelated/family * autosomal/X-chromosomal variants
// need careful tests before any change
class MetaScoreTest : public ModelFitter {
 public:
  MetaScoreTest() : model(NULL), modelAuto(NULL), modelX(NULL), useBolt(false) {
    this->modelName = "MetaScore";
    af = -1.;
    fitOK = false;
    useFamilyModel = false;
    isHemiRegion = false;
    headerOutputted = false;
    indexResult = true;
    outputSE = false;
  }
  virtual ~MetaScoreTest() {
    if (modelAuto) {
      delete modelAuto;
      modelAuto = NULL;
    }
    if (modelX) {
      delete modelX;
      modelX = NULL;
    }
  }
  virtual int setParameter(const ModelParser& parser) {
    this->outputGwama = parser.hasTag("gwama");
    this->outputSE = parser.hasTag("se");
    this->useBolt = parser.hasTag("bolt");
    return 0;
  }
  // fitting model
  virtual int fit(DataConsolidator* dc) {
    Matrix& genotype = dc->getGenotype();
    return this->fitWithGivenGenotype(genotype, dc);
  }
  int fitWithGivenGenotype(Matrix& genotype, DataConsolidator* dc) {
    this->useFamilyModel = dc->hasKinship();

    // check column name for hemi region
    this->isHemiRegion = dc->isHemiRegion(0);

    // fit null model for the header
    if (isHemiRegion) {
      model = modelX;
      if (!model) {
        model = modelX = createModel(this->useFamilyModel, isBinaryOutcome());
        model->hemiRegion = true;
      }
    } else {
      model = modelAuto;
      if (!model) {
        model = modelAuto =
            createModel(this->useFamilyModel, isBinaryOutcome());
        model->hemiRegion = false;
      }
    }
    if (!model) return -1;

    // calculate site-based statise
    if (isHemiRegion) {
      // use female only
      dc->countRawGenotypeFromFemale(0, &counter);
    } else {
      dc->countRawGenotype(0, &counter);
    }
    af = counter.getAF();

    // for binary trait, also count by case and controls
    if (isBinaryOutcome()) {
      if (isHemiRegion) {
        // use female only
        dc->countRawGenotypeFromFemaleCase(0, &caseCounter);
        dc->countRawGenotypeFromFemaleControl(0, &ctrlCounter);
      } else {
        dc->countRawGenotypeFromCase(0, &caseCounter);
        dc->countRawGenotypeFromControl(0, &ctrlCounter);
      }
    }
    // place this after calculate site statistics e.g. AF
    // fit null model
    if (model->needToFitNullModel || dc->isPhenotypeUpdated() ||
        dc->isCovariateUpdated()) {
      // copyCovariateAndIntercept(genotype.rows, covariate, &cov);
      fitOK = (0 == model->FitNullModel(genotype, dc));
      if (!fitOK) return -1;
      model->needToFitNullModel = false;
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
    fitOK = (0 == model->TestCovariate(genotype, dc));

    if (useFamilyModel) {
      this->af = model->GetAF(dc->getGenotype(), dc);
    }
    return (fitOK ? 0 : -1);
  }

  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    //  In the header part, successfully estimated null model will be
    //  printed. So we have to defer printing header in the fit() function.
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
    if (outputSE) {
      result.addHeader("ALT_EFFSIZE_SE");
    }
    result.addHeader("PVALUE");
    return;
  }

  void writeSummaryAndHeader(FileWriter* fp, const Result& siteInfo) {
    if (g_SummaryHeader) {
      g_SummaryHeader->outputHeader(fp);
    }
    if (model) {
      model->PrintNullModel(fp, g_SummaryHeader->getCovLabel());
    }
    if (outputGwama) {
      // output
      model->DoExtraWork(fp, siteInfo, this);
    }

    siteInfo.writeHeaderTab(fp);
    result.writeHeaderLine(fp);
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    if (!headerOutputted) {
      writeSummaryAndHeader(fp, siteInfo);
      headerOutputted = true;
    }

    siteInfo.writeValueTab(fp);

    result.clearValue();
    if (af >= 0.0) {
      if (!isBinaryOutcome()) {
        result.updateValue("AF", af);
      } else {
        static char afString[128];
        snprintf(afString, 128, "%g:%g:%g", af, caseCounter.getAF(),
                 ctrlCounter.getAF());
        result.updateValue("AF", afString);
      }
    }

    if (!isBinaryOutcome()) {
      result.updateValue("INFORMATIVE_ALT_AC", counter.getAC());
      result.updateValue("CALL_RATE", counter.getCallRate());
      result.updateValue("HWE_PVALUE", counter.getHWE());
      result.updateValue("N_REF", counter.getNumHomRef());
      result.updateValue("N_HET", counter.getNumHet());
      result.updateValue("N_ALT", counter.getNumHomAlt());
    } else {
      static const int buffLen = 128;
      static char buff[buffLen];

      snprintf(buff, 128, "%g:%g:%g", counter.getAC(), caseCounter.getAC(),
               ctrlCounter.getAC());
      result.updateValue("INFORMATIVE_ALT_AC", buff);

      snprintf(buff, 128, "%g:%g:%g", counter.getCallRate(),
               caseCounter.getCallRate(), ctrlCounter.getCallRate());
      result.updateValue("CALL_RATE", buff);

      snprintf(buff, 128, "%g:%g:%g", counter.getHWE(), caseCounter.getHWE(),
               ctrlCounter.getHWE());
      result.updateValue("HWE_PVALUE", buff);

      snprintf(buff, 128, "%d:%d:%d", counter.getNumHomRef(),
               caseCounter.getNumHomRef(), ctrlCounter.getNumHomRef());
      result.updateValue("N_REF", buff);
      snprintf(buff, 128, "%d:%d:%d", counter.getNumHet(),
               caseCounter.getNumHet(), ctrlCounter.getNumHet());
      result.updateValue("N_HET", buff);
      snprintf(buff, 128, "%d:%d:%d", counter.getNumHomAlt(),
               caseCounter.getNumHomAlt(), ctrlCounter.getNumHomAlt());
      result.updateValue("N_ALT", buff);
    }

    if (fitOK) {
      const double u = model->GetU();
      const double v = model->GetV();
      result.updateValue("U_STAT", u);
      result.updateValue("SQRT_V_STAT", sqrt(v));
      result.updateValue("ALT_EFFSIZE", model->GetEffect());
      if (outputSE && v > 0.) {
        result.updateValue("ALT_EFFSIZE_SE", model->GetEffectSE());
      }
      result.updateValue("PVALUE", model->GetPvalue());
    }
    result.writeValueLine(fp);
  }

 private:
  class MetaBase {
   public:
    MetaBase() : needToFitNullModel(true) {}
    virtual ~MetaBase() {}
    virtual int FitNullModel(Matrix& genotype, DataConsolidator* dc) = 0;
    virtual int TestCovariate(Matrix& genotype, DataConsolidator* dc) = 0;
    virtual double GetAF(Matrix& geno, DataConsolidator* dc) = 0;
    virtual void PrintNullModel(FileWriter* fp,
                                const std::vector<std::string>& covLabel) = 0;
    // this method is for printing extra stuffs
    virtual int DoExtraWork(FileWriter* fp, const Result& siteInfo,
                            ModelFitter* model) {
      return 0;
    };
    virtual double GetU() = 0;
    virtual double GetV() = 0;
    virtual double GetEffect() = 0;
    virtual double GetEffectSE() = 0;
    virtual double GetPvalue() = 0;

    bool needToFitNullModel;
    bool hemiRegion;

   protected:
    Matrix cov;
  };  // class MetaBase
  class MetaFamQtl : public MetaBase {
   public:
    MetaFamQtl() : lmm(FastLMM::SCORE, FastLMM::MLE) {}
    int FitNullModel(Matrix& genotype, DataConsolidator* dc) {
      Matrix& phenotype = dc->getPhenotype();
      Matrix& covariate = dc->getCovariate();
      copyCovariateAndIntercept(genotype.rows, covariate, &cov);

      bool fitOK;
      if (!hemiRegion) {
        fitOK = (0 ==
                 lmm.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(),
                                  *dc->getKinshipSForAuto()));
      } else {
        if (!dc->hasKinshipForX()) {
          fitOK = false;
          return -1;
        }
        fitOK = (0 ==
                 lmm.FitNullModel(cov, phenotype, *dc->getKinshipUForX(),
                                  *dc->getKinshipSForX()));
      }
      if (fitOK) {
        needToFitNullModel = false;
        return 0;
      }
      return -1;
    }
    int TestCovariate(Matrix& genotype, DataConsolidator* dc) {
      Matrix& phenotype = dc->getPhenotype();
      Matrix& cov = dc->getCovariate();
      if (!hemiRegion) {
        return lmm.TestCovariate(cov, phenotype, genotype,
                                 *dc->getKinshipUForAuto(),
                                 *dc->getKinshipSForAuto());
      } else {
        return lmm.TestCovariate(cov, phenotype, genotype,
                                 *dc->getKinshipUForX(),
                                 *dc->getKinshipSForX());
      }
    }
    double GetAF(Matrix& geno, DataConsolidator* dc) {
      if (!hemiRegion) {
        return lmm.FastGetAF(*dc->getKinshipUForAuto(),
                             *dc->getKinshipSForAuto(), dc->getGenotype());
      } else {
        return lmm.FastGetAF(*dc->getKinshipUForX(), *dc->getKinshipSForX(),
                             dc->getGenotype());
      }
    }
    void PrintNullModel(FileWriter* fp,
                        const std::vector<std::string>& covLabel) {
      Vector beta;
      lmm.GetNullCovEst(&beta);
      Matrix betaSd;
      lmm.GetNullCovB(&betaSd);
      double sigmaG2 = lmm.GetSigmaG2();
      double sigmaE2 = lmm.GetSigmaE2();

      fp->write("##NullModelEstimates\n");
      fp->write("## - Name\tBeta\tSD\n");
      fp->printf("## - Intercept\t%g\t%g\n", beta[0], betaSd[0][0]);
      const int n = covLabel.size();
      for (int i = 0; i < n; ++i) {
        if (i + 1 >= beta.Length()) break;
        fp->printf("## - %s\t%g\t%g\n", covLabel[i].c_str(), beta[i + 1],
                   betaSd[i + 1][i + 1]);
      }
      // sigma
      fp->printf("## - SigmaG2\t%g\tNA\n", sigmaG2);
      fp->printf("## - SigmaE2\t%g\tNA\n", sigmaE2);
    }
    int DoExtraWork(FileWriter* fp, const Result& siteInfo,
                    ModelFitter* model) {
      std::string s = model->getPrefix();
      s += ".MetaScore.assoc.factors";
      FileWriter fout(s.c_str());
      double sigmaK = lmm.GetSigmaK();
      double sigma1 = lmm.GetSigma1();
      fout.printf("Sigma_K\t%g\n", sigmaK);
      fout.printf("Sigma_1\t%g\n", sigma1);

      s = model->getPrefix();
      s += ".MetaScore.assoc.zy";
      Matrix m;
      lmm.GetCovZY(&m);
      FileWriter fw(s.c_str());
      for (int i = 0; i < m.rows; ++i) {
        for (int j = 0; j < m.cols; ++j) {
          if (j) fw.write("\t");
          fw.printf("%g", m[i][j]);
        }
        fw.write("\n");
      }
      return 0;
    }
    double GetU() { return lmm.GetUStat(); }
    double GetV() { return lmm.GetVStat(); }
    double GetEffect() {
      return lmm.GetVStat() != 0.0 ? lmm.GetUStat() / lmm.GetVStat() : 0.0;
    }
    double GetEffectSE() { return lmm.GetSE(); }
    double GetPvalue() { return lmm.GetPvalue(); }

   private:
    FastLMM lmm;
  };  // class MetaFamQtl
  class MetaUnrelatedQtl : public MetaBase {
   public:
    int FitNullModel(Matrix& genotype, DataConsolidator* dc) {
      Matrix& phenotype = dc->getPhenotype();
      Matrix& covariate = dc->getCovariate();

      copyCovariateAndIntercept(genotype.rows, covariate, &this->cov);
      copyPhenotype(phenotype, &this->pheno);

      bool fitOK = linear.FitNullModel(cov, pheno);
      if (!fitOK) return -1;
      needToFitNullModel = false;
      sigma2 = linear.GetSigma2();
      return 0;
    }
    int TestCovariate(Matrix& genotype, DataConsolidator* dc) {
      bool fitOK = linear.TestCovariate(cov, pheno, genotype);
      if (!fitOK) return -1;
      return 0;
    }
    double GetAF(Matrix& geno, DataConsolidator* dc) {
      assert(false);
      ;  // should not reach here
      return 0.0;
    }
    void PrintNullModel(FileWriter* fp,
                        const std::vector<std::string>& covLabel) {
      Vector& beta = linear.GetNullCovEst();
      Matrix& betaSd = linear.GetNullCovB();

      fp->write("##NullModelEstimates\n");
      fp->write("## - Name\tBeta\tSD\n");
      fp->printf("## - Intercept\t%g\t%g\n", beta[0], betaSd[0][0]);
      const int n = covLabel.size();
      for (int i = 0; i < n; ++i) {
        if (i + 1 >= beta.Length()) break;
        fp->printf("## - %s\t%g\t%g\n", covLabel[i].c_str(), beta[i + 1],
                   betaSd[i + 1][i + 1]);
      }
      // sigma
      fp->printf("## - Sigma2\t%g\tNA\n", sigma2);
    }
    double GetU() { return linear.GetU()[0][0] / sigma2; }
    double GetV() { return linear.GetV()[0][0] / sigma2 / sigma2; }
    double GetEffect() {
      return linear.GetV()[0][0] != 0.0 ? linear.GetBeta()[0][0] : 0.0;
    }
    double GetEffectSE() { return linear.GetSEBeta(0); }
    double GetPvalue() { return linear.GetPvalue(); }

   private:
    LinearRegressionScoreTest linear;
    double sigma2;
    Vector pheno;
  };  //   class MetaUnrelatedQtl
  class MetaFamBinary : public MetaBase {
   public:
    MetaFamBinary() : lmm(FastLMM::SCORE, FastLMM::MLE) {
      lmm.disableCenterGenotype();
    }
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
                 lmm.FitNullModel(cov, phenotype, *dc->getKinshipUForAuto(),
                                  *dc->getKinshipSForAuto()));
      } else {
        if (!dc->hasKinshipForX()) {
          fitOK = false;
          return -1;
        }
        fitOK = (0 ==
                 lmm.FitNullModel(cov, phenotype, *dc->getKinshipUForX(),
                                  *dc->getKinshipSForX()));
      }
      if (fitOK) {
        needToFitNullModel = false;
        return 0;
      }
      return -1;
    }
    int TestCovariate(Matrix& genotype, DataConsolidator* dc) {
      Matrix& phenotype = dc->getPhenotype();
      Matrix& cov = dc->getCovariate();
      if (!hemiRegion) {
        return lmm.TestCovariate(cov, phenotype, genotype,
                                 *dc->getKinshipUForAuto(),
                                 *dc->getKinshipSForAuto());
      } else {
        return lmm.TestCovariate(cov, phenotype, genotype,
                                 *dc->getKinshipUForX(),
                                 *dc->getKinshipSForX());
      }
    }
    double GetAF(Matrix& geno, DataConsolidator* dc) {
      if (!hemiRegion) {
        return lmm.FastGetAF(*dc->getKinshipUForAuto(),
                             *dc->getKinshipSForAuto(), dc->getGenotype());
      } else {
        return lmm.FastGetAF(*dc->getKinshipUForX(), *dc->getKinshipSForX(),
                             dc->getGenotype());
      }
    }
    void PrintNullModel(FileWriter* fp,
                        const std::vector<std::string>& covLabel) {
      Vector beta;
      lmm.GetNullCovEst(&beta);
      Matrix betaSd;
      lmm.GetNullCovB(&betaSd);
      double sigmaG2 = lmm.GetSigmaG2();
      double sigmaE2 = lmm.GetSigmaE2();

      fp->write("##NullModelEstimates\n");
      fp->write("## - Name\tBeta\tSD\n");
      fp->printf("## - Intercept\t%g\t%g\n", beta[0], betaSd[0][0]);
      const int n = covLabel.size();
      for (int i = 0; i < n; ++i) {
        if (i + 1 >= beta.Length()) break;
        fp->printf("## - %s\t%g\t%g\n", covLabel[i].c_str(), beta[i + 1],
                   betaSd[i + 1][i + 1]);
      }
      // sigma
      fp->printf("## - SigmaG2\t%g\tNA\n", sigmaG2);
      fp->printf("## - SigmaE2\t%g\tNA\n", sigmaE2);
    }
    double GetU() { return lmm.GetUStat() * this->b; }
    double GetV() { return lmm.GetVStat() * this->b * this->b; }
    double GetEffect() {
      if (lmm.GetVStat() != 0.0 && this->b != 0.0) {
        return lmm.GetUStat() / lmm.GetVStat() / b;
      }
      return 0.0;
    }
    double GetEffectSE() {
      const double v = GetV();
      return v != 0.0 ? 1.0 / sqrt(v) / b : 0.0;
    }
    double GetPvalue() { return lmm.GetPvalue(); }

   private:
    void calculateB();

   private:
    FastLMM lmm;
    double alpha;
    double b;
  };  // class MetaFamBinary
  class MetaUnrelatedBinary : public MetaBase {
   public:
    MetaUnrelatedBinary() : useMLE(false) {}
    int FitNullModel(Matrix& genotype, DataConsolidator* dc) {
      Matrix& phenotype = dc->getPhenotype();
      Matrix& covariate = dc->getCovariate();

      copyCovariateAndIntercept(genotype.rows, covariate, &this->cov);
      copyPhenotype(phenotype, &this->pheno);

      if (!useMLE) {
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
        /*
          no need to adjust for logistic regression
        if (nCtrl > 0) {
          alpha = log(1.0 * nCase / nCtrl);
        } else {
          alpha = 500.;
        }
        calculateB();
        */
      }
      // fit null model
      bool fitOK = logistic.FitNullModel(cov, pheno, 100);
      if (!fitOK) return -1;
      needToFitNullModel = false;
      return 0;
    }
    int TestCovariate(Matrix& genotype, DataConsolidator* dc) {
      bool fitOK = logistic.TestCovariate(cov, pheno, genotype);
      if (!fitOK) return -1;

      if (useMLE) {
        // this part may be optimized by using approximations
        // fit alternative model
        Matrix& covariate = dc->getCovariate();
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
      return 0;
    }
    double GetAF(Matrix& geno, DataConsolidator* dc) {
      assert(false);
      ;  // should not reach here
      return 0.0;
    }
    void PrintNullModel(FileWriter* fp,
                        const std::vector<std::string>& covLabel) {
      Vector& beta = logistic.GetNullCovEst();
      Matrix& betaSd = logistic.GetNullCovB();

      fp->write("##NullModelEstimates\n");
      fp->write("## - Name\tBeta\tSD\n");
      fp->printf("## - Intercept\t%g\t%g\n", beta[0], betaSd[0][0]);
      const int n = covLabel.size();
      for (int i = 0; i < n; ++i) {
        if (i + 1 >= beta.Length()) break;
        fp->printf("## - %s\t%g\t%g\n", covLabel[i].c_str(), beta[i + 1],
                   betaSd[i + 1][i + 1]);
      }
      // sigma
      fp->printf("## - Sigma2\tNA\tNA\n");
    }
    double GetU() { return logistic.GetU()[0][0]; }
    double GetV() { return logistic.GetV()[0][0]; }
    double GetEffect() {
      if (!useMLE) {
        if (logistic.GetU()[0][0] != 0.0) {
          return logistic.GetU()[0][0] / logistic.GetV()[0][0];
        }
      } else {
        return logisticAlt.GetCovEst()[1];
      }
      return 0.0;
    }
    double GetEffectSE() {
      const double v = GetV();
      return v != 0.0 ? 1.0 / sqrt(v) : 0.0;
    }
    double GetPvalue() { return logistic.GetPvalue(); }

    // private:
    //  void calculateB();

   private:
    LogisticRegressionScoreTest logistic;
    Vector pheno;
    Matrix X;  // intercept, cov(optional) and genotype
    // double alpha;
    // double b;

    bool useMLE;
    LogisticRegression logisticAlt;
  };  // class MetaUnrelatedBinary

  class MetaFamQtlBolt : public MetaBase {
   public:
    MetaFamQtlBolt() { fprintf(stderr, "MetaFamQtlBolt model started\n"); }
    int FitNullModel(Matrix& genotype, DataConsolidator* dc) {
      const std::string& fn = dc->getBoltGenotypeFilePrefix();

      // fit null model
      bool fitOK = 0 == bolt_.FitNullModel(fn, &dc->getPhenotype());
      if (!fitOK) return -1;
      needToFitNullModel = false;
      return 0;
    }
    int TestCovariate(Matrix& genotype, DataConsolidator* dc) {
      bool fitOK = 0 == bolt_.TestCovariate(genotype);
      if (!fitOK) return -1;
      return 0;
    }
    double GetAF(Matrix& geno, DataConsolidator* dc) { return bolt_.GetAF(); }
    void PrintNullModel(FileWriter* fp,
                        const std::vector<std::string>& covLabel) {}
    double GetU() { return bolt_.GetU(); }
    double GetV() { return bolt_.GetV(); }
    double GetEffect() { return bolt_.GetEffect(); }
    double GetEffectSE() {
      const double v = bolt_.GetV();
      ;
      return v != 0.0 ? 1.0 / sqrt(v) : 0.0;
    }
    double GetPvalue() { return bolt_.GetPvalue(); }

   private:
    BoltLMM bolt_;
  };  // class MetaFamQtlBolt
  MetaBase* createModel(bool familyModel, bool binaryOutcome) {
    MetaBase* ret = NULL;
    if (this->useBolt) {
      if (binaryOutcome) {
        fprintf(stderr, "BoltLMM does not support binary outcomes! Exit...\n");
        exit(1);
      }
      ret = new MetaFamQtlBolt;
      return ret;
    }

    if (familyModel && !binaryOutcome) {
      ret = new MetaFamQtl;
    }
    if (familyModel && binaryOutcome) {
      ret = new MetaFamBinary;
    }
    if (!familyModel && binaryOutcome) {
      ret = new MetaUnrelatedBinary;
    }
    if (!familyModel && !binaryOutcome) {
      ret = new MetaUnrelatedQtl;
    }
    return ret;
  }

  void reset() {
    ModelFitter::reset();
    counter.reset();
    caseCounter.reset();
    ctrlCounter.reset();
  }

 private:
  MetaBase* model;
  MetaBase* modelAuto;
  MetaBase* modelX;

 protected:
  // allow inherited-class to change
  bool useBolt;

 private:
  double af;  // overall af (unadjust or adjusted by family structure)

  bool fitOK;

  GenotypeCounter counter;
  GenotypeCounter caseCounter;
  GenotypeCounter ctrlCounter;

  bool useFamilyModel;
  bool isHemiRegion;  // is the variant tested in hemi region?
  bool headerOutputted;
  bool outputGwama;
  bool outputSE;
};  // MetaScoreTest

class MetaDominantTest : public MetaScoreTest {
 public:
  MetaDominantTest() : MetaScoreTest() { this->modelName = "MetaDominant"; }
  virtual int fit(DataConsolidator* dc) {
    dc->codeGenotypeForDominantModel(&geno);
    return fitWithGivenGenotype(geno, dc);
  }

 private:
  Matrix geno;
};

class MetaRecessiveTest : public MetaScoreTest {
 public:
  MetaRecessiveTest() : MetaScoreTest() { this->modelName = "MetaRecessive"; }
  virtual int fit(DataConsolidator* dc) {
    dc->codeGenotypeForRecessiveModel(&geno);
    return fitWithGivenGenotype(geno, dc);
  }

 private:
  Matrix geno;
};

class MetaScoreBoltTest : public MetaScoreTest {
 public:
  MetaScoreBoltTest() : MetaScoreTest() { this->modelName = "MetaScoreBolt"; }
  virtual int setParameter(const ModelParser& parser) {
    MetaScoreTest::setParameter(parser);
    this->useBolt = true;
    return 0;
  }
};

class MetaCovBase;
class MetaCovTest : public ModelFitter {
 public:
  typedef int Genotype;
  typedef int Covariate;
  struct Pos {
    std::string chrom;
    int pos;
  };
  struct Loci {
    Pos pos;
    Genotype geno;
    Covariate covXZ;
    // Genotype geno;
    // std::vector<float> covXZ;  // cov(geno, covariate)
  };

 public:
  MetaCovTest(int windowSize);
  virtual ~MetaCovTest();
  virtual int setParameter(const ModelParser& parser) {
    this->outputGwama = parser.hasTag("gwama");
    return 0;
  }
  // fitting model
  virtual int fit(DataConsolidator* dc) {
    Matrix& genotype = dc->getGenotype();
    return this->fitWithGivenGenotype(genotype, dc);
  }
  int fitWithGivenGenotype(Matrix& genotype, DataConsolidator* dc);

  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    if (g_SummaryHeader) {
      g_SummaryHeader->outputHeader(fp);
    }

    result.writeHeaderLine(fp);
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    this->fout = fp;
    while (queue.size() && getWindowSize(queue, loci) > windowSize) {
      printCovariance(fout, queue, isBinaryOutcome());
      genoPool.deallocate(queue.front().geno);
      genoCovPool.deallocate(queue.front().covXZ);
      queue.pop_front();
    }
    if (fitOK) {
      queue.push_back(loci);
      ++numVariant;
    }
    // result.writeValueLine(fp);
  }

 private:
  void assignGenotype(Matrix& genotype, Genotype& genoIdx);
  /**
   * @return max integer if different chromosome; or return difference between
   * head and tail locus.
   */
  int getWindowSize(const std::deque<Loci>& loci, const Loci& newOne) {
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
  }
  // X is symmetric matrix
  // return: a' * X * b
  // = \sum_i \sum_j a_i * X_{ij} * b_j
  // = \sum_i {a_i * X_ii * b_i + \sum_{j!=i} (a_i * b_j + a_j *b_i) * X_{ij}
  float computeQuadraticForm(FloatMatrixRef& a, Matrix& X, FloatMatrixRef& b) {
    // const int n = X.rows;
    // float s = 0.;
    // for (int i = 0; i < n; ++i) {
    //   s += a[i] * X[i][i] * b[i];
    //   for (int j = 0; j < i; ++j) {
    //     s += (a[i] * b[j] + a[j] * b[i]) * X[i][j];
    //   }
    // }
    // return s;

    Eigen::MatrixXf XE;
    G_to_Eigen(X, &XE);
    REF_TO_EIGEN(a, aE);
    REF_TO_EIGEN(b, bE);
    return (aE.transpose() * XE * bE)(0, 0);
  }
  // return = covX1X2 - covX1Z' * covZZInv * covX2Z
  float computeScaledXX(const float covX1X2, FloatMatrixRef& covX1Z,
                        FloatMatrixRef& covX2Z, Matrix& covZZInv) {
    float ret = covX1X2;
    if (covX1Z.ncol_) {
      ret -= computeQuadraticForm(covX1Z, covZZInv, covX2Z);
    }
    // fprintf(stderr, "%d: ret = %g\n", __LINE__, ret);
    return ret;
  }
  /**
   * @return 0
   * print the covariance for the front of loci to the rest of loci
   */
  int printCovariance(FileWriter* fp, const std::deque<Loci>& lociQueue,
                      bool binaryOutcome);
  void appendToString(const std::vector<int>& position, std::string* out) {
    std::string& s = *out;
    for (size_t i = 0; i < position.size(); ++i) {
      if (i) s += ',';
      s += toString(position[i]);
    }
  }
  void appendToString(const std::vector<float>& vec, const float scale,
                      std::string* out) {
    std::string& s = *out;
    for (size_t i = 0; i < vec.size(); ++i) {
      if (i) s += ',';
      s += toString(vec[i] * scale);
    }
  }
  void appendToString(FloatMatrixRef& vec, const float scale,
                      std::string* out) {
    const int n = vec.nrow_;
    std::string& s = *out;
    for (int i = 0; i < n; ++i) {
      if (i) s += ',';
      s += toString(vec.memory_[i] * scale);
    }
  }
  void appendToString(Matrix& mat, const float scale, std::string* out) {
    if (mat.cols != mat.rows) {
      fprintf(stderr, "only square matrix is supported!\n");
    }
    std::string& s = *out;
    for (int i = 0; i < mat.rows; ++i) {
      for (int j = 0; j <= i; ++j) {
        if (i || j) s += ',';
        s += floatToString(mat[i][j] * scale);
      }
    }
  }

 private:
  MetaCovBase* createModel(bool familyModel, bool binaryOutcome);

 private:
  MetaCovBase* model;
  MetaCovBase* modelAuto;
  MetaCovBase* modelX;

 protected:
  // allow inherited-class to change
  bool useBolt;

 private:
  std::deque<Loci> queue;
  RingMemoryPool genoPool;     // store genotypes
  RingMemoryPool genoCovPool;  // store G'Z , e.g. genotype * covariate)
  int numVariant;
  int nSample;
  int nCovariate;
  // Vector mleVarY;  // variance term of Y_i for i = 1,..., N th sample
  FileWriter* fout;
  int windowSize;
  Loci loci;
  bool fitOK;
  std::vector<float> covXX;
  Matrix covZZ;
  Matrix covZZInv;
  bool useFamilyModel;
  Matrix cov;
  bool isHemiRegion;  // is the variant tested in hemi region
  std::vector<int> position;
  bool outputGwama;
};  // MetaCovTest

class MetaDominantCovTest : public MetaCovTest {
 public:
  MetaDominantCovTest(int windowSize) : MetaCovTest(windowSize) {
    this->modelName = "MetaDominantCov";
  }
  virtual int fit(DataConsolidator* dc) {
    dc->codeGenotypeForDominantModel(&geno);
    return fitWithGivenGenotype(geno, dc);
  }

 private:
  Matrix geno;
};

class MetaRecessiveCovTest : public MetaCovTest {
 public:
  MetaRecessiveCovTest(int windowSize) : MetaCovTest(windowSize) {
    this->modelName = "MetaRecessiveCov";
  }
  virtual int fit(DataConsolidator* dc) {
    dc->codeGenotypeForRecessiveModel(&geno);
    return fitWithGivenGenotype(geno, dc);
  }

 private:
  Matrix geno;
};

class MetaCovBoltTest : public MetaCovTest {
 public:
  MetaCovBoltTest(int windowSize) : MetaCovTest(windowSize) {
    this->modelName = "MetaCovBolt";
  }
  virtual int setParameter(const ModelParser& parser) {
    MetaCovTest::setParameter(parser);
    MetaCovTest::useBolt = true;
    return 0;
  }
};

#if 0
class MetaSkewTest: public ModelFitter{
 private:
  typedef std::vector<double> Genotype;
  struct Pos{
    std::string chrom;
    int pos;
  }
    struct Loci{
      Pos pos;
      Genotype geno;
    }

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
  }
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
        fprintf(stderr, "Sample size changed at [ %s:%s ]\n", siteInfo["CHROM"].c_str(), siteInfo["POS"].c_str());
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
    }
    loci.geno.resize(nSample);
    for (int i = 0; i < nSample; ++i) {
      loci.geno[i] = genotype[i][0];
    }
    fitOK = true;
    return (fitOK ? 0 : -1);
  }

  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    if (g_SummaryHeader) {
      g_SummaryHeader->outputHeader(fp);
    }

    /* siteInfo.writeHeaderTab(fp); */
    // fprintf(fp, "AF\tStat\tDirection\tPvalue\n");
    result.writeHeaderLine(fp);
  }
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
    }
    if (fitOK) {
      queue.push_back(loci);
      ++numVariant;
    }
    // result.writeValueLine(fp);
  }

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
  }
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
  }

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
  }
    struct Loci{
      Pos pos;
      Genotype geno;
    }

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
  }
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
        fprintf(stderr, "Sample size changed at [ %s:%s ]\n", siteInfo["CHROM"].c_str(), siteInfo["POS"].c_str());
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
    }
    loci.geno.resize(nSample);
    for (int i = 0; i < nSample; ++i) {
      loci.geno[i] = genotype[i][0];
    }
    fitOK = true;
    return (fitOK ? 0 : -1);
  }

  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    if (g_SummaryHeader) {
      g_SummaryHeader->outputHeader(fp);
    }

    /* siteInfo.writeHeaderTab(fp); */
    // fprintf(fp, "AF\tStat\tDirection\tPvalue\n");
    result.writeHeaderLine(fp);
  }
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
    }
    if (fitOK) {
      queue.push_back(loci);
      ++numVariant;
    }
    // result.writeValueLine(fp);
  }

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
  }
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
  }

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

#define MULTIPLE_TRAIT_SCORE_TEST_BLOCK_SIZE 8
class MultipleTraitScoreTest : public ModelFitter {
 public:
  MultipleTraitScoreTest()
      : nSample(-1),
        linear(MULTIPLE_TRAIT_SCORE_TEST_BLOCK_SIZE),
        fitOK(false),
        needToFitNullModel(true),
        numResult(0),
        blockSize(MULTIPLE_TRAIT_SCORE_TEST_BLOCK_SIZE),
        sites(MULTIPLE_TRAIT_SCORE_TEST_BLOCK_SIZE) {
    this->modelName = "MultipleTraitScore";
  }
  ~MultipleTraitScoreTest() { flushOutput(); }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();
    const FormulaVector& tests = *dc->getFormula();

    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }
    nSample = genotype.rows;

    // af = getMarkerFrequency(genotype, 0);
    // if (isMonomorphicMarker(genotype, 0)) {
    //   fitOK = false;
    //   return -1;
    // }
    // copyCovariateAndIntercept(genotype.rows, covariate, &cov);
    // copyPhenotype(phenotype, &this->pheno);

    if (!isBinaryOutcome()) {  // qtl
      if (needToFitNullModel || dc->isPhenotypeUpdated() ||
          dc->isCovariateUpdated()) {
        fitOK = linear.FitNullModel(covariate, phenotype, tests);
        if (!fitOK) {
          warnOnce("Multiple trait score test failed when fitting null model.");
          return -1;
        }
        needToFitNullModel = false;
      }
      fitOK = linear.AddGenotype(genotype);
    } else {
      warnOnce(
          "Multiple trait score test model does not support binary trait yet.");
      return -1;
      // if (needToFitNullModel || dc->isPhenotypeUpdated() ||
      //     dc->isCovariateUpdated()) {
      //   fitOK = logistic.FitNullModel(cov, pheno, 100);
      //   if (!fitOK) {
      //     warnOnce("Single variant score test failed in fitting null
      //     model.");
      //     return -1;
      //   }
      //   calculateConstant(phenotype);
      //   needToFitNullModel = false;
      // }
      // fitOK = logistic.TestCovariate(cov, pheno, genotype);
    }
    return (fitOK ? 0 : -1);
  }

  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    // result.addHeader("AF");
    // result.addHeader("U");
    // result.addHeader("V");
    // result.addHeader("STAT");
    // result.addHeader("DIRECTION");
    // result.addHeader("EFFECT");
    // result.addHeader("SE");
    result.addHeader("U_STAT");
    result.addHeader("V_STAT");
    result.addHeader("PVALUE");
    result.writeHeaderLine(fp);
    this->fp = fp;
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    // siteInfo.writeValueTab(fp);
    // result.clearValue();
    // result.updateValue("AF", af);
    sites[numResult].clear();
    siteInfo.writeValueTab(&sites[numResult]);
    numResult++;
    if (fitOK) {
      if (!isBinaryOutcome()) {
        // formatValue(linear.GetU(), &ustat);
        // formatValue(linear.GetV(), &vstat);
        // formatValue(linear.GetPvalue(), &pvalue);

        // result.updateValue("U_STAT", ustat);
        // result.updateValue("V_STAT", vstat);
        // result.updateValue("PVALUE", pvalue);
        if (numResult == blockSize) {
          flushOutput();
        }
        // const double u = linear.GetU()[0][0];
        // const double v = linear.GetV()[0][0];
        // result.updateValue("U", u);
        // result.updateValue("V", v);
        // result.updateValue("STAT", linear.GetStat());
        // if (u != 0) {
        //   result.updateValue("DIRECTION", linear.GetU()[0][0] > 0 ? "+" :
        //   "-");
        // }
        // if (v > 0) {
        //   result.updateValue("EFFECT", linear.GetBeta()[0][0]);
        //   result.updateValue("SE", 1.0 / sqrt(v));
        // }
        // result.updateValue("PVALUE", linear.GetPvalue());
      } else {
        // const double u = logistic.GetU()[0][0];
        // const double v = logistic.GetV()[0][0];
        // result.updateValue("U", u);
        // result.updateValue("V", v);
        // result.updateValue("STAT", logistic.GetStat());
        // if (u != 0) {
        //   result.updateValue("DIRECTION",
        //                      logistic.GetU()[0][0] > 0 ? "+" : "-");
        // }
        // if (v > 0 && b > 0) {
        //   result.updateValue("EFFECT", u / v / b);
        //   result.updateValue("SE", 1.0 / sqrt(v) / b);  // need to verify
        // }
        // result.updateValue("PVALUE", logistic.GetPvalue());
      }
    }
    // result.writeValueLine(fp);
  }
  void flushOutput() {
    fitOK = linear.TestCovariateBlock();
    if (fitOK) {
      if (!isBinaryOutcome()) {
        for (int i = 0; i < numResult; ++i) {
          fp->write(sites[i]);

          formatValue(linear.GetU(i), &ustat);
          formatValue(linear.GetV(i), &vstat);
          formatValue(linear.GetPvalue(i), &pvalue);

          result.updateValue("U_STAT", ustat);
          result.updateValue("V_STAT", vstat);
          result.updateValue("PVALUE", pvalue);

          result.writeValueLine(fp);
        }
      }
    }
    linear.flush();
    numResult = 0;
  }
  void formatValue(const Vector& v, std::string* out) {
    const int n = v.Length();
    (*out).clear();
    for (int i = 0; i < n; ++i) {
      if (i) {
        (*out) += ",";
      }
      (*out) += toString(v[i]);
    }
  }

 private:
  // double b;  // a constant
  // double af;
  int nSample;
  // Vector pheno;
  MultipleTraitLinearRegressionScoreTest linear;
  bool fitOK;
  bool needToFitNullModel;
  int numResult;
  int blockSize;
  std::string ustat;
  std::string vstat;
  std::string pvalue;
  FileWriter* fp;
  std::vector<std::string> sites;
  // Matrix cov;
};  // MultipleTraitScoreTest

class FastMultipleTraitScoreTest : public ModelFitter {
 public:
  FastMultipleTraitScoreTest()
      : nSample(-1),
        linear(MULTIPLE_TRAIT_SCORE_TEST_BLOCK_SIZE),
        fitOK(false),
        needToFitNullModel(true),
        numResult(0),
        blockSize(MULTIPLE_TRAIT_SCORE_TEST_BLOCK_SIZE),
        sites(MULTIPLE_TRAIT_SCORE_TEST_BLOCK_SIZE) {
    this->modelName = "FastMultipleTraitScore";
  }
  ~FastMultipleTraitScoreTest() { flushOutput(); }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();
    const FormulaVector& tests = *dc->getFormula();

    if (genotype.cols != 1) {
      fitOK = false;
      return -1;
    }
    nSample = genotype.rows;

    // af = getMarkerFrequency(genotype, 0);
    // if (isMonomorphicMarker(genotype, 0)) {
    //   fitOK = false;
    //   return -1;
    // }
    // copyCovariateAndIntercept(genotype.rows, covariate, &cov);
    // copyPhenotype(phenotype, &this->pheno);

    if (!isBinaryOutcome()) {  // qtl
      if (needToFitNullModel || dc->isPhenotypeUpdated() ||
          dc->isCovariateUpdated()) {
        fitOK = linear.FitNullModel(covariate, phenotype, tests);
        if (!fitOK) {
          warnOnce("Multiple trait score test failed when fitting null model.");
          return -1;
        }
        needToFitNullModel = false;
      }
      fitOK = linear.AddGenotype(genotype);
    } else {
      warnOnce(
          "Multiple trait score test model does not support binary trait yet.");
      return -1;
      // if (needToFitNullModel || dc->isPhenotypeUpdated() ||
      //     dc->isCovariateUpdated()) {
      //   fitOK = logistic.FitNullModel(cov, pheno, 100);
      //   if (!fitOK) {
      //     warnOnce("Single variant score test failed in fitting null
      //     model.");
      //     return -1;
      //   }
      //   calculateConstant(phenotype);
      //   needToFitNullModel = false;
      // }
      // fitOK = logistic.TestCovariate(cov, pheno, genotype);
    }
    return (fitOK ? 0 : -1);
  }

  // write result header
  void writeHeader(FileWriter* fp, const Result& siteInfo) {
    siteInfo.writeHeaderTab(fp);
    // result.addHeader("AF");
    // result.addHeader("U");
    // result.addHeader("V");
    // result.addHeader("STAT");
    // result.addHeader("DIRECTION");
    // result.addHeader("EFFECT");
    // result.addHeader("SE");
    result.addHeader("U_STAT");
    result.addHeader("V_STAT");
    result.addHeader("PVALUE");
    result.writeHeaderLine(fp);
    this->fp = fp;
  }
  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    // siteInfo.writeValueTab(fp);
    // result.clearValue();
    // result.updateValue("AF", af);
    sites[numResult].clear();
    siteInfo.writeValueTab(&sites[numResult]);
    numResult++;
    if (fitOK) {
      if (!isBinaryOutcome()) {
        // formatValue(linear.GetU(), &ustat);
        // formatValue(linear.GetV(), &vstat);
        // formatValue(linear.GetPvalue(), &pvalue);

        // result.updateValue("U_STAT", ustat);
        // result.updateValue("V_STAT", vstat);
        // result.updateValue("PVALUE", pvalue);
        if (numResult == blockSize) {
          flushOutput();
        }
        // const double u = linear.GetU()[0][0];
        // const double v = linear.GetV()[0][0];
        // result.updateValue("U", u);
        // result.updateValue("V", v);
        // result.updateValue("STAT", linear.GetStat());
        // if (u != 0) {
        //   result.updateValue("DIRECTION", linear.GetU()[0][0] > 0 ? "+" :
        //   "-");
        // }
        // if (v > 0) {
        //   result.updateValue("EFFECT", linear.GetBeta()[0][0]);
        //   result.updateValue("SE", 1.0 / sqrt(v));
        // }
        // result.updateValue("PVALUE", linear.GetPvalue());
      } else {
        // const double u = logistic.GetU()[0][0];
        // const double v = logistic.GetV()[0][0];
        // result.updateValue("U", u);
        // result.updateValue("V", v);
        // result.updateValue("STAT", logistic.GetStat());
        // if (u != 0) {
        //   result.updateValue("DIRECTION",
        //                      logistic.GetU()[0][0] > 0 ? "+" : "-");
        // }
        // if (v > 0 && b > 0) {
        //   result.updateValue("EFFECT", u / v / b);
        //   result.updateValue("SE", 1.0 / sqrt(v) / b);  // need to verify
        // }
        // result.updateValue("PVALUE", logistic.GetPvalue());
      }
    }
    // result.writeValueLine(fp);
  }
  void flushOutput() {
    fitOK = linear.TestCovariateBlock();
    if (fitOK) {
      if (!isBinaryOutcome()) {
        for (int i = 0; i < numResult; ++i) {
          fp->write(sites[i]);

          formatValue(linear.GetU(i), &ustat);
          formatValue(linear.GetV(i), &vstat);
          formatValue(linear.GetPvalue(i), &pvalue);

          result.updateValue("U_STAT", ustat);
          result.updateValue("V_STAT", vstat);
          result.updateValue("PVALUE", pvalue);

          result.writeValueLine(fp);
        }
      }
    }
    linear.flush();
    numResult = 0;
  }
  void formatValue(const Vector& v, std::string* out) {
    const int n = v.Length();
    (*out).clear();
    for (int i = 0; i < n; ++i) {
      if (i) {
        (*out) += ",";
      }
      (*out) += toString(v[i]);
    }
  }

 private:
  // double b;  // a constant
  // double af;
  int nSample;
  // Vector pheno;
  FastMultipleTraitLinearRegressionScoreTest linear;
  bool fitOK;
  bool needToFitNullModel;
  int numResult;
  int blockSize;
  std::string ustat;
  std::string vstat;
  std::string pvalue;
  FileWriter* fp;
  std::vector<std::string> sites;
  // Matrix cov;
};  // FastMultipleTraitScoreTest

class DumpModel : public ModelFitter {
 public:
  DumpModel(const char* prefix) {
    this->prefix = prefix;
    this->modelName = "DumpData";
  }
  // write result header
  void writeHeader(FileWriter* fp,
                   const Result& siteInfo) {  // e.g. column headers.
    siteInfo.writeHeaderTab(fp);
    fp->write("FileName\n");

    this->header = siteInfo.joinHeader();
  }
  // fitting model
  int fit(DataConsolidator* dc) {
    Matrix& phenotype = dc->getPhenotype();
    Matrix& genotype = dc->getGenotype();
    Matrix& covariate = dc->getCovariate();

    this->phenotype = phenotype;
    this->genotype = genotype;
    this->covariate = covariate;
    // copy genotype column label
    for (int i = 0; i < genotype.cols; ++i) {
      this->genotype.SetColumnLabel(i, genotype.GetColumnLabel(i));
    }
    // copy covaraite column label
    for (int i = 0; i < covariate.cols; ++i) {
      this->covariate.SetColumnLabel(i, covariate.GetColumnLabel(i));
    }
    // copy row labels
    this->rowLabel = dc->getRowLabel();
    return 0;
  }

  // write model output
  void writeOutput(FileWriter* fp, const Result& siteInfo) {
    std::string fn = this->prefix + "\t" + siteInfo.joinValue() + ".data";

    for (size_t i = 0; i < fn.size(); ++i) {
      if (fn[i] == '\t') fn[i] = '.';
    }

    siteInfo.writeValueTab(fp);
    fp->printf("%s\n", fn.c_str());

    // write header
    FILE* fDump = fopen(fn.c_str(), "wt");
    fprintf(fDump, "ID\t");
    fprintf(fDump, "%s\t", this->header.c_str());
    if (phenotype.cols == 1) {
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
    }
    for (int i = 0; i < covariate.cols; i++) {
      if (strlen(covariate.GetColumnLabel(i)) == 0) {
        fprintf(fDump, "\tC%d", i);
      } else {
        fprintf(fDump, "\t%s", covariate.GetColumnLabel(i));
      }
    }
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
    }
    fclose(fDump);
  }

  void reset() {}  // for particular class to call when fitting repeatedly

 private:
  Matrix phenotype;
  Matrix genotype;
  Matrix covariate;
  std::string prefix;
  std::string header;
  std::vector<std::string> rowLabel;
};  // end DumpModel

#endif /* _MODELFITTER_H_ */
