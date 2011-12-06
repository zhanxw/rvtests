#ifndef _MODELFITTER_H_
#define _MODELFITTER_H_

#include "VCFData.h"
#include "regression/LogisticRegression.h"
#include "regression/LinearRegression.h"

// take X, Y, Cov and fit model
class ModelFitter{
public:
    // write result header
    virtual void writeHeader(FILE* fp) = 0;
    // fitting model
    virtual int fit(Matrix* geno,
                    Vector* pheno,
                    Matrix* cov, 
                    FILE* fp) = 0;
}; // end ModelFitter

class LogisticModelFitter: public ModelFitter{
  public:
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "PVALUE\n");
    };
    // fitting model
    int fit(Matrix* geno,
            Vector* pheno,
            Matrix* cov, 
            FILE* fp) {
        LogisticRegressionScoreTest lrst;
        bool fitOK = false;
        if (cov) {
            fitOK = lrst.FitNullModel(*cov, *pheno, 100);
        }
        fitOK = lrst.TestCovariate(*cov, *pheno, *geno);
        fprintf(fp, "%lf\n", lrst.getPvalue());
    };
}; // LogisticModelFitter
#if 0
class LinearModelFitter: public ModelFitter{
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "PVALUE\n");
    };
    // fitting model
    int fit(Matrix* geno,
            Vector* pheno,
            Matrix* cov, 
            FILE* fp) {
        LinearRegressionScoreTest lrst;
        lrst.FitLogisticModel( *cov, *pheno, *geno, 100);
        fprintf(fp, "%lf\n", lr.getPvalue());
    };
};
#endif
#endif /* _MODELFITTER_H_ */
