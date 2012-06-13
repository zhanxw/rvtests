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

// take X, Y, Cov and fit model
// note, ModelFitter will use VCFData as READ-ONLY data structure,
// and collapsing results are stored internally.
class ModelFitter{
public:
  // write result header
  virtual void writeHeader(FILE* fp, const char* prependString) = 0;
  // fitting model
  virtual int fit(Matrix& phenotype, Matrix& genotype) = 0;
  // write model output
  virtual void writeOutput(FILE* fp, const char* prependString) = 0;

  ModelFitter(){
    this->modelName = "Unassigned_Model_Name";
    this->binaryOutcome = false; // default: using continuous outcome
  };
  const std::string& getModelName() const { return this->modelName; };
  virtual void reset() {}; // for particular class to call when fitting repeatedly
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
};

class SingleVariantWaldTest: public ModelFitter{
public:
  SingleVariantWaldTest(){
    this->modelName = "SingleWald";
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fprintf(fp, "Beta\tSE\tPvalue\n");
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype) {
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    copyPhenotype(phenotype, &this->Y);
    copyGenotypeWithIntercept(genotype, &this->X);

    if (!isBinaryOutcome()) {
      fitOK = linear.FitLinearModel(this->X, this->Y);
    } else {
      fitOK = logistic.FitLogisticModel(this->X, this->Y, 100);
    }
    return (fitOK ? 0 : 1);
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (!fitOK) {
      fputs("NA\tNA\tNA\n", fp);
    } else {
      if (!isBinaryOutcome()) {
        double se = sqrt(linear.GetCovB()[1][1]);
        fprintf(fp, "%.3lf\t%.3lf\t%.3lf\n", linear.GetCovEst()[1], se, linear.GetAsyPvalue()[1]);
      } else {
        double se = sqrt(logistic.GetCovB()[0][0]);
        fprintf(fp, "%.3lf\t%.3lf\t%.3lf\n", logistic.GetCovEst()[1], se, logistic.GetAsyPvalue()[1]);
      }
    }
  };
private:
  Matrix X; // 1 + geno
  Vector Y; // phenotype
  LinearRegression linear;
  LogisticRegression logistic;
  int nCov;
  int nGeno;
  bool fitOK;
}; // SingleVariantWaldTest

class SingleVariantScoreTest: public ModelFitter{
public:
  SingleVariantScoreTest(){
    this->modelName = "SingleScore";
  };

  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fprintf(fp, "NSample\tAF\tPvalue\n");
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype) {
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    nSample = genotype.rows;
    af = getMarkerFrequency(genotype, 0);
    Matrix intercept;
    intercept.Dimension(genotype.rows, 1);
    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i< phenotype.rows; i++){
      intercept[i][0] = 1.0;
      pheno[i] = phenotype[i][0];
    }
    fitOK = lrst.FitNullModel(intercept, pheno);
    if (!fitOK) return -1;
    fitOK = lrst.TestCovariate(intercept, pheno, genotype);
    return (fitOK ? 0 : -1);
  };
#if 0
  int fit(VCFData& data) {
    if (data.covariate && data.covariate->cols > 0) {
      fitOK = lrst.FitNullModel( *data.covariate, *data.extractPhenotype(), 100);
      if (!fitOK)
        return -1;

      fitOK = lrst.TestCovariate( *data.covariate, *data.extractPhenotype(), *data.collapsedGenotype);
      if (!fitOK)
        return -1;
    } else {
      fitOK = lrst.TestCovariate( *data.collapsedGenotype, *data.extractPhenotype());
      if (!fitOK)
        return -1;
    }
    return 0;
  };
#endif
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (fitOK)
      fprintf(fp, "%d\t%.3f\t%g\n", nSample, af, lrst.GetPvalue());
    else
      fputs("NA\n", fp);
  };
private:
  double af;
  int nSample;
  LinearRegressionScoreTest lrst;
  bool fitOK;
}; // SingleVariantScoreTest

class SingleVariantLogisticWaldTest: public ModelFitter{
public:
  SingleVariantLogisticWaldTest(){
    this->modelName = "SingleWald";
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fprintf(fp, "Beta\tSE\tPvalue");
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype) {
    cbind(&X, &genotype, NULL, true);
    fitOK = lr.FitLogisticModel(X, phenotype, 100);
    return (fitOK ? 0 : 1);
  };

#if 0
  int fit(VCFData& data) {
    if (data.covariate && data.covariate->cols > 0)
      cbind(&X, data.genotype, data.covariate, true);
    else
      cbind(&X, data.genotype, NULL, true);

    fitOK = lr.FitLogisticModel(X, *data.extractPhenotype(), 100);
    return (fitOK ? 0 : 1);
  };
#endif
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (fitOK) {
      double se = sqrt(lr.GetCovB()[0][0]);
      fprintf(fp, "%.3lf\t%.3lf\t%.3lf", lr.GetCovEst()[0], se, lr.GetAsyPvalue()[0]);
    } else{
      fputs("NA\tNA\tNA", fp);
    }
  };
private:
  void cbind(Matrix* out, Matrix* a, Matrix* b, bool addIntercept) {
    assert(out && a);

    int totalCol = a->cols + (b ? b->cols : 0) + (addIntercept ? 1 : 0);
    if (b)
      assert ( a->rows == b->rows);
    out->Dimension(a->rows, totalCol);

    for (int r = 0; r < a->rows; r ++ ){
      int beginCol = 0;
      for (int c = 0; c < a->cols; c++) {
        (*out)[r][c] = (*a)[r][c];
      }
      beginCol += a->cols;

      if (b) {
        for (int c = 0; c < b->cols; c++) {
          (*out)[r][c] = (*a)[r][beginCol + c];
        }
        beginCol += b->cols;
      }

      if (addIntercept) {
        (*out)[r][beginCol] = 1.0;
      }
    }
  };
private:
  Matrix X; // geno + 1 + covariate
  LogisticRegression lr;
  bool fitOK;

}; // SingleVariantLogisticWaldTest

class SingleVariantLogisticScoreTest: public ModelFitter{
public:
  SingleVariantLogisticScoreTest(){
    this->modelName = "SingleScore";
  };

  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fprintf(fp, "Pvalue\n");
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype) {
    Matrix intercept;
    intercept.Dimension(genotype.rows, 1);
    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i< phenotype.rows; i++){
      intercept[i][0] = 1.0;
      pheno[i] = phenotype[i][0];
    }
    fitOK = lrst.FitNullModel(intercept, pheno, 100);
    if (!fitOK) return -1;
    fitOK = lrst.TestCovariate(intercept, pheno, genotype);
    return (fitOK ? 0 : -1);
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (fitOK)
      fprintf(fp, "%.3f\n", lrst.GetPvalue());
    else
      fputs("NA\n", fp);
  };
private:
  LogisticRegressionScoreTest lrst;
  bool fitOK;
}; // SingleVariantLogisticScoreTest


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
  int fit(Matrix& phenotype, Matrix& genotype) {
    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    Matrix intercept;
    intercept.Dimension(genotype.rows, 1);
    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i< phenotype.rows; i++){
      intercept[i][0] = 1.0;
      pheno[i] = phenotype[i][0];
    }

    cmcCollapse(genotype, &collapsedGenotype);

    if (isBinaryOutcome()) {
      fitOK = logistic.FitNullModel(intercept, pheno, 100);
      if (!fitOK) return -1;
      fitOK = logistic.TestCovariate(intercept, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    } else {
      fitOK = linear.FitNullModel(intercept, pheno);
      if (!fitOK) return -1;
      fitOK = linear.TestCovariate(intercept, pheno, collapsedGenotype);
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
  int fit(Matrix& phenotype, Matrix& genotype) {
    if (genotype.cols == 0) {
      fitOK = false;
      return 1;
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
  int fit(Matrix& phenotype, Matrix& genotype) {
    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    Matrix intercept;
    intercept.Dimension(genotype.rows, 1);
    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i< phenotype.rows; i++){
      intercept[i][0] = 1.0;
      pheno[i] = phenotype[i][0];
    }

    zegginiCollapse(genotype, &collapsedGenotype);

    if (isBinaryOutcome()) {
      fitOK = logistic.FitNullModel(intercept, pheno, 100);
      if (!fitOK) return -1;
      fitOK = logistic.TestCovariate(intercept, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    } else {
      fitOK = linear.FitNullModel(intercept, pheno);
      if (!fitOK) return -1;
      fitOK = linear.TestCovariate(intercept, pheno, collapsedGenotype);
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
MadsonBrowningTest(int nPerm): nPerm(nPerm) {
    this->modelName = "MadsonBrowningTest";
  }
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (isBinaryOutcome())
      fprintf(fp, "NumVariant\tNumPerm\tActualPerm\tStat\tNumX\tPval\n");
    else
      fprintf(fp, "NumVariant\tNA\n");
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype) {
    if (!isBinaryOutcome()) {
      fitOK = false;
      return -1;
    }

    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    Matrix intercept;
    intercept.Dimension(genotype.rows, 1);
    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i< phenotype.rows; i++){
      intercept[i][0] = 1.0;
      pheno[i] = phenotype[i][0];
    }

    madsonBrowningCollapse(genotype, pheno, &collapsedGenotype);


    fitOK = logistic.FitNullModel(intercept, pheno, 100);
    if (!fitOK) return -1;
    fitOK = logistic.TestCovariate(intercept, pheno, collapsedGenotype);
    if (!fitOK) return -1;
    // record observed stat
    this->stat = logistic.GetStat(); // a chi-dist

    // permutation part
    int failed = 0;
    for (int i = 0; i < this->nPerm; ++i) {
      permute(&pheno);
      fitOK = logistic.TestCovariate(intercept, pheno, collapsedGenotype);
      if (!fitOK) {
        if (failed < 10) {
          failed++;
          continue;
        }else {
          fitOK = false;
          return -1;
        }
      }
      ++ this->actualPerm;
      // record new stats
      double pStat = logistic.GetStat();
      if (pStat > this->stat) {
        this->numX ++;
      }
    } // end permutation
    return (fitOK ? 0 : -1);
  };
  void reset() {
    this->actualPerm = 0;
    this->numX = 0;
  }
  double getPvalue() {
    return (1.0 + this->numX) / (1.0 + this->actualPerm);
  }
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (isBinaryOutcome()) {
      if (fitOK){
        fprintf(fp, "%d\t%d\t%d\t%g\t%d\t%g\n",
                this->numVariant,
                this->nPerm,
                this->actualPerm,
                this->stat,
                this->numX,
                this->getPvalue());
      } else {
        fprintf(fp, "%d\tNA\tNA\tNA\tNA\tNA\n", this->numVariant);
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
  int nPerm;
  int numX;
  int actualPerm;
  double stat;
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
  int fit(Matrix& phenotype, Matrix& genotype) {
    this->numVariant = genotype.cols;
    if (genotype.cols == 0) {
      fitOK = false;
      return -1;
    }
    Matrix intercept;
    intercept.Dimension(genotype.rows, 1);
    Vector pheno;
    pheno.Dimension(phenotype.rows);
    for (int i = 0; i< phenotype.rows; i++){
      intercept[i][0] = 1.0;
      pheno[i] = phenotype[i][0];
    }

    fpCollapse(genotype, &collapsedGenotype);

    if (isBinaryOutcome()) {
      fitOK = logistic.FitNullModel(intercept, pheno, 100);
      if (!fitOK) return -1;
      fitOK = logistic.TestCovariate(intercept, pheno, collapsedGenotype);
      return (fitOK ? 0 : -1);
    } else {
      fitOK = linear.FitNullModel(intercept, pheno);
      if (!fitOK) return -1;
      fitOK = linear.TestCovariate(intercept, pheno, collapsedGenotype);
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

/**
 * Implementation of Alkes Price's VT
 */
class VariableThresholdPrice: public ModelFitter{
public:
VariableThresholdPrice(int nPerm): nPerm(nPerm) {
    this->modelName = "VariableThresholdPrice";
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype) {
    rearrangeGenotypeByFrequency(genotype, &sortedGenotype, &this->freq);
    convertToReferenceAlleleCount(&sortedGenotype);
    collapseGenotype(&sortedGenotype);
    transpose(&sortedGenotype); // now each row is a collapsed genoype at certain frequency cutoff
    copyPhenotype(phenotype, &this->phenotype);

    this->zmax = 100.0;
    double z;

    if (isBinaryOutcome()) {
      for (int i = 0; i < sortedGenotype.rows; ++i) {
        z = calculateZForBinaryTrait(this->phenotype, this->sortedGenotype[i]);
        if ( z > this->zmax) {
          zmax = z;
          this->optimalFreq = freq[i];
        }
      }

      // begin permutation
      for (int i = 0; i < nPerm; ++i) {
        permute(&this->phenotype);
        ++actualPerm;
        for (int j = 0; j < sortedGenotype.rows; ++j){
          z = calculateZForBinaryTrait(this->phenotype, this->sortedGenotype[j]);
          if (z > this->zmax) {
            ++ this->numX;
          };
        }
      };

    } else {
      centerVector(&this->phenotype);
      for (int i = 0; i < sortedGenotype.rows; ++i) {
        z = calculateZForContinuousTrait(this->phenotype, this->sortedGenotype[i]);
        if ( z > this->zmax){
          zmax = z;
          this->optimalFreq = freq[i];
        }
      }

      // begin permutation
      for (int i = 0; i < nPerm; ++i) {
        permute(&this->phenotype);
        ++actualPerm;
        for (int j = 0; j < sortedGenotype.rows; ++j){
          z = calculateZForContinuousTrait(this->phenotype, this->sortedGenotype[j]);
          if (z > this->zmax) {
            ++ this->numX;
          };
        }
      };
    };

    fitOK = true;
    return 0;
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fputs("NumPerm\tActualPerm\tOptimalFreq\tZmax\tNumX\tPvalue\n", fp);
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fprintf(fp, "%d\t%d\t%g\t%g\t%d\t%g\n",
            this->nPerm,
            this->actualPerm,
            this->optimalFreq,
            this->zmax,
            this->numX,
            this->calculatePvalue()
            );
  };
  void reset() {
    fitOK = true;
    this->zmax = 100.0;
    this->pvalue = 100.0;
    this->actualPerm = 0;
    this->numX = 0;
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
  double calculatePvalue() {
    return (1.0 + numX) / (1.0 + actualPerm);
  };
private:
  Matrix sortedGenotype;
  std::vector<double> freq;
  bool fitOK;
  Vector phenotype;
  double zmax;
  double pvalue;
  int nPerm;
  int actualPerm;
  int numX; // x is the number of permutation for which Zmax is higher in unpermuted data than permutated data.
  double optimalFreq; // the frequency cutoff in unpermutated data which give smallest pvalue
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
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    model[0].writeHeader(fp, "FreqThreshold\t");
  };
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype) {

    if (genotype.cols > modelLen) {
      resize(genotype.cols);
      reset();
    }

    rearrangeGenotypeByFrequency(genotype, &sortedGenotype, &this->freq);
    for (int i = genotype.cols - 1; i >=0; --i){
      sortedGenotype.Dimension( genotype.rows, i + 1);
      if ( model[i].fit(phenotype, sortedGenotype) ) {
        fitOK = false;
      }
    }

    return 0;
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
      fprintf(fp, "%.3lf", outFreq[i]);
    }
    fputc('\t', fp);
    for (int i = 0; i < pvalue.size(); i++) {
      if (i) fputc(',', fp);
      if (pvalue[i] > 0)
        fprintf(fp, "%.3lf", pvalue[i]);
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
  SkatTest() {
    this->modelName = "Skat";
  };
  // write result header
  void writeHeader(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    fprintf(fp, "NMARKER\tQ\tPvalue\n");
  };
  void reset() {
    this->skat.Reset();
  }
  // fitting model
  int fit(Matrix& phenotype, Matrix& genotype) {
    // fill it wegith
    // NOTE: our frequency calculation is slightly different than SKAT, so we will need to adjust it back
    weight.Dimension(genotype.cols);
    for (int i = 0; i < weight.Length(); i++) {
      double freq = getMarkerFrequency(genotype, i);
      freq = (freq * (genotype.rows * 2 + 1) -1 ) / (genotype.rows *2 );
      if (freq > 1e-30) { // avoid dividing zero
        weight[i] =  gsl_ran_beta_pdf(freq, 1.0, 25.0);  /// use beta(MAF, 1, 25)
        weight[i] *= weight[i];
        // fprintf(stderr, "weight(%d, %d, %f ) = %f\n", 1, 25, freq, weight[i]);
      } else {
        weight[i] = 0.0;
      }
    };

    // get ynulll
    // ynull is mean of y (removing genotypes) in model Ynull ~ X (aka Ynull ~ X + 0.0 * G )
    X.Dimension(genotype.rows, 1);
    for (int i = 0; i < genotype.rows; ++i) {
      X[i][0] = 1.0;
    }
    if (isBinaryOutcome()) {
      fitOK = logistic.FitLogisticModel(X, phenotype, 100);
      if (!fitOK) {
        return -1;
      }
      ynull = logistic.GetPredicted();
      v = logistic.GetVariance();
    } else {
      fitOK = linear.FitLinearModel(X, phenotype);
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
    skat.CalculatePValue(phenotype, ynull, X, v, genotype, weight);
    this->pValue = skat.GetPvalue();
    return 0;
  };
  // write model output
  void writeOutput(FILE* fp, const char* prependString) {
    fputs(prependString, fp);
    if (!fitOK){
      fprintf(fp, "%d\tNA\tNA\n", this->weight.Length());
    } else {
      if (isBinaryOutcome() ) {
        fprintf(fp, "%d\t%g\t%g\n", this->weight.Length(), this->skat.GetQ(), this->pValue);
      } else {
        fprintf(fp, "%d\t%g\t%g\n", this->weight.Length(), this->skat.GetQ(), this->pValue);
      }
    }
  };

private:
  Matrix X; // n by (p+1) matrix, people by covariate (note intercept is needed);
  Vector v;
  Vector weight;
  LogisticRegression logistic;
  LinearRegression linear;
  Vector ynull;
  Skat skat;
  bool fitOK;
  double pValue;
  int nMarker;
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
  int fit(Matrix& phenotype, Matrix& genotype) {
    if (!isBinaryOutcome()) {
      fitOK = false;
      return -1;
    }

    this->resize(genotype.rows, genotype.cols);
    this->nn = this->nPerm;
    this->qq = 1;
    this->aa = 0.0;
    this->mafUpper = 1.0;

    for (int j = 0; j < genotype.cols; ++j) {
      for (int i = 0; i < genotype.rows; ++i) {
        xdatIn[j * genotype.cols + i] = genotype[i][j];
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
    do_kbac_test(&this->pValue, &twosided);

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
    if (numPeople > this->ylen) {
      delete[] this->ydatIn;
      this->ydatIn = new double[ numPeople];
      this->ylen = numPeople;
      resized = true;
    }
    if (numMarker > this->xcol) {
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
  int fit(Matrix& phenotype, Matrix& genotype) {
    this->phenotype = phenotype;
    this->genotype = genotype;
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
      fprintf(fDump, "\n");
    };
    fclose(fDump);

  };

  void reset() {
  }; // for particular class to call when fitting repeatedly
private:
  Matrix phenotype;
  Matrix genotype;
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
  double freq = 1.0 * (ac + 1) / (an + 1);
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
