#ifndef _MODELFITTER_H_
#define _MODELFITTER_H_

#include "VCFData.h"
#include "regression/LogisticRegression.h"
#include "regression/LinearRegression.h"

// take X, Y, Cov and fit model
class ModelFitter{
public:
    // write result header
    virtual void outputHeader(FILE* fp) = 0;
    // fitting model
    virtual int fit(VCFData& data, Collapsor& c, FILE* fp) = 0;

    ModelFitter(){
        this->modelName = "Unassigned Model Name";
    };
    const std::string& getModelName() { return this->modelName; };
    void reset() {}; // for particular class to call
protected:
    std::string modelName;
}; // end ModelFitter

class CaseControlFreqSummary: public ModelFitter{
  public:
    CaseControlFreqSummary(){
        this->modelName = "FrequencySummary";
    };
    void outputHeader(FILE* fp){
        fprintf(fp, "NumpTotalSample\tNumSNP\tNumCaseTotal\tNumCase0\tNumCase1\tNumCase2\tNumCaseMissing\tNumControlTotal\tNumControl0\tNumControl1\tNumControl2\tNumControlMissing");
    };
    int fit(VCFData& data, Collapsor& c, FILE* fp){
        int freq[2][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}}; // case
        const int CASE = 0;
        const int CONTROL = 1;
        int group = -1;
        for (int p = 0; p < data.genotype->cols; p++){
            int pheno = (int) ((*data.phenotype)[p][0]);
            switch (pheno){
            case 0:
                group = CONTROL;
                break;
            case 1:
                group = CASE;
                break;
            default:
                //treat as missing
                continue;
            }
            for (int m = 0; m < data.genotype->rows; m++){
                int geno = (int) ((*data.genotype)[m][p]);
                switch (geno){
                case 0:
                case 1:
                case 2:
                    freq[group][geno] ++ ;
                    break;
                case -9:
                    freq[group][3] ++;
                    break;
                default:
                    fprintf(stderr, "unrecognized genotypes: %.2f\n", (*data.genotype)[m][p]);
                    break;
                }
            }
        }
        int nCaseTotal = freq[CASE][0] + freq[CASE][1] + freq[CASE][2] + freq[CASE][3];
        int nControlTotal = freq[CONTROL][0] + freq[CONTROL][1] + freq[CONTROL][2] + freq[CONTROL][3];
        fprintf(fp, "%d\t%d", data.genotype->cols, data.genotype->rows);
        fprintf(fp, "\t%d\t%d\t%d\t%d", nCaseTotal, freq[CASE][0], freq[CASE][1], freq[CASE][2], freq[CASE][3]);
        fprintf(fp, "\t%d\t%d\t%d\t%d", nControlTotal, freq[CONTROL][0], freq[CONTROL][1], freq[CONTROL][2], freq[CONTROL][3]);
    };
}; // end class CaseControlFreqSummary

class LogisticModelScoreTest: public ModelFitter{
public:
    LogisticModelScoreTest(){
        this->modelName = "LogisticModelScoreTest";
        this->dirtyMark = true;
    };
    // write result header
    void outputHeader(FILE* fp) {
        fprintf(fp, "%s.PVALUE", this->modelName.c_str());
    };
    // fitting model
    int fit(VCFData& d, Collapsor& c, FILE* fp) {
        c.collapseGenotype(&d, NULL);
        Matrix* geno = d.collapsedGenotype;
        Vector* pheno = d.extractPhenotype();
        Matrix* cov = d.covariate;

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

class VariableThresholdTest: public ModelFitter{
  public:
    VariableThresholdTest(){
        this->modelName = "VT";
    };

    // write result header
    void outputHeader(FILE* fp){
        fprintf(fp, "ThresholdCutoffs\tPvalues");
    };
    // fitting model
    int fit(VCFData& d, Collapsor& c, FILE* fp){
        int numPeople = d.people2Idx.size();
        int numMarker = d.marker2Idx.size();
        
        d.calculateFrequency(VCFData::FREQ_CONTORL_ONLY);
        std::map<double, int> freqCutOff;
        std::map<double, int>::const_iterator freqCutOff_iter;
        for (int m = 0; m < numMarker; m++){
            if (d.markerFreq[m] > 1e-6) {
                freqCutOff[d.markerFreq[m]] = m;
            }
        }
        std::vector<double> pvalue;

        // collapsor.setCollapsingStrategy(Collapsor::MADSON_BROWNING);
        this->progressiveCollapsor.setCollapsingStrategy(Collapsor::PROGRESSIVE);
        int param[2]; 
        param[1] = Collapsor::CMC;
        for (freqCutOff_iter = freqCutOff.begin();
             freqCutOff_iter != freqCutOff.end();
             freqCutOff_iter ++){

            param[0] = freqCutOff_iter->second;
            this->progressiveCollapsor.collapseGenotype(&d, &param);
            pt.FitLinearModel(*d.collapsedGenotype, *d.extractPhenotype(), 1000, 0.10); // 0.10 = 0.05 * 2 ( try use adaptive permutation here).
            pvalue.push_back(pt.getPvalue());
        }
        // output results;
        for (freqCutOff_iter = freqCutOff.begin();
             freqCutOff_iter != freqCutOff.end();
             freqCutOff_iter ++){
            if (freqCutOff_iter != freqCutOff.begin()) fprintf(fp, ",");
            fprintf(fp, "%.3f", freqCutOff_iter->first);
        }
        fprintf(fp, "\t");
        for (int i = 0; i < pvalue.size(); i++){
            if (i) fprintf(fp, ",");
            fprintf(fp, "%.3f", pvalue[i]);
        }
    };
  private:
    Collapsor progressiveCollapsor;
    LinearPermutationTest pt;
};




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
