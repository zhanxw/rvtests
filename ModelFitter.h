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
    const std::string& getModelName() { return this->modelName; };
protected:
    std::string modelName;
}; // end ModelFitter

class LogisticModelScoreTest: public ModelFitter{
public:
    LogisticModelScoreTest(){
        this->modelName = "LogisticModelScoreTest";
        this->dirtyMark = true;
    };
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "%s.PVALUE", this->modelName.c_str());
    };
    // fitting model
    int fit(Matrix* geno,
            Vector* pheno,
            Matrix* cov,
            FILE* fp) {

        bool fitOK = false;
        // when cov exists, we may cache null model
        if (cov && cov->cols > 0) {


            // deal with dirtyMark
            if (this->phenoCache != pheno ||
                this->covCache != cov)
                this->dirtyMark = true;

            if (dirtyMark) {
                fitOK = lrst.FitNullModel(*cov, *pheno, 100);
                this->dirtyMark = false;
                this->phenoCache = pheno;
                this->covCache = cov;
            }

            if (!fitOK) {
                fputs("-1.00", fp);
                return -1;
            }

            fitOK = lrst.TestCovariate(*cov, *pheno, *geno);
            if (!fitOK) {
                fputs("-1.00", fp);
                return -1;
            }

            fprintf(fp, "%lf\n", lrst.getPvalue());
            return 0;
        } else {
            // no covariate
            fitOK = lrst.TestCovariate(*geno, *pheno);
            if (!fitOK) {
                fputs("-1.00", fp);
                return -1;
            }
            fprintf(fp, "%lf\n", lrst.getPvalue());
            return 0;
        }
        //  should not reach here.
        return -2;
    };
private:
    LogisticRegressionScoreTest lrst;
    bool dirtyMark; // dirtyMark == true: need to refit Null Model
    void* phenoCache;
    void* covCache;
}; // LogisticModelScoreTest

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
