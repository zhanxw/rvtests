#ifndef _MODELFITTER_H_
#define _MODELFITTER_H_

#include "VCFData.h"
#include "regression/LogisticRegression.h"
#include "regression/LinearRegression.h"

// take X, Y, Cov and fit model
// note, ModelFitter will use VCFData as READ-ONLY data structure, 
// and collapsing results are stored internally.

class ModelFitter{
public:
    // write result header
    virtual void writeHeader(FILE* fp) = 0;
    // fitting model
    virtual int fit(VCFData& data) = 0;
    // write model output
    virtual void writeOutput(FILE* fp) = 0;

    ModelFitter(){
        this->modelName = "Unassigned_Model_Name";
    };
    const std::string& getModelName() const { return this->modelName; };
    void reset() {}; // for particular class to call
protected:
    std::string modelName;
}; // end ModelFitter

class SingleVariantHeader: public ModelFitter{
  public:
    SingleVariantHeader() {
        modelName = "SingleVariantModelHeader";
        this->reset();
    };
    void writeHeader(FILE* fp) {
        fprintf(fp, "NumpTotalSample\tNumSNP\tNumCaseTotal\tNumCase0\tNumCase1\tNumCase2\tNumCaseMissing\tNumControlTotal\tNumControl0\tNumControl1\tNumControl2\tNumControlMissing");
    };
    int fit (VCFData& data) {
        this->numPeople = data.genotype->cols;
        this->numMarker = data.genotype->rows;
        int group = -1;
        for (int p = 0; p < numPeople; p++){
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
        nCaseTotal = freq[CASE][0] + freq[CASE][1] + freq[CASE][2] + freq[CASE][3];
        nControlTotal = freq[CONTROL][0] + freq[CONTROL][1] + freq[CONTROL][2] + freq[CONTROL][3];
    };
    void writeOutput(FILE* fp) {
        fprintf(fp, "%d\t%d", numPeople, numMarker);
        fprintf(fp, "\t%d\t%d\t%d\t%d", nCaseTotal, freq[CASE][0], freq[CASE][1], freq[CASE][2], freq[CASE][3]);
        fprintf(fp, "\t%d\t%d\t%d\t%d", nControlTotal, freq[CONTROL][0], freq[CONTROL][1], freq[CONTROL][2], freq[CONTROL][3]);
    };
    void reset() {
        freq[0][0] = freq[0][1] = freq[0][2] = freq[0][3] = 0; // case 
        freq[1][0] = freq[1][1] = freq[1][2] = freq[1][3] = 0; // control
        nCaseTotal = nControlTotal = 0;
    };
  private:
    const static int CASE = 0;
    const static int CONTROL = 1;

    int numPeople;
    int numMarker;
    int freq[2][4];
    int nCaseTotal;
    int nControlTotal;
}; //SingleVariantHeader

class CollapsingHeader: public ModelFitter{
    public:
    // write result header
    void writeHeader(FILE* fp) {
        fputs("NumVariant", fp);
    };
    // fitting model
    int fit(VCFData& data) {
        this->numVariant = data.collapsedGenotype.cols;
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
        fprintf(fp, "%d", this->numVariant);
    };
    private:
    int numVariant = 0;
}; // CollapsingHeader

class SingleVariantWaldTest: public ModelFitter{
    public:
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "Beta\tSE\tPvalue");
    };
    // fitting model
    int fit(VCFData& data) {
        lr.FitModel();
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
    };
    private:
    Matrix X; // geno + 1 + covariate
    LogisticRegression lr;
}; // SingleVariantWaldTest

class SingleVariantScoreTest: public ModelFitter{
    public:
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "Pvalue");
    };
    // fitting model
    int fit(VCFData& data) {
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
        fprintf(fp, "%.3f", lrst.GetPvalue());
    };
    private:
    LogisticRegressionScoreTest lrst;
}; // SingleVariantScoreTest

class CMCTest: public ModelFitter{
    public:
    // write result header
    void writeHeader(FILE* fp) {
    };
    // fitting model
    int fit(VCFData& data) {
        // collapsing
        // fit model
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
    };
    prviate:
    Matrix g;
    LogisticRegressionScoreTest lrst;
}; // CMCTest

class ZegginiTest: public ModelFitter{
    // write result header
    void writeHeader(FILE* fp) {
    };
    // fitting model
    int fit(VCFData& data) {
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
    };
}; // ZegginiTest

class MadsonBrowningTest: public ModelFitter{
    // write result header
    void writeHeader(FILE* fp) {
    };
    // fitting model
    int fit(VCFData& data) {
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
    };
}; // MadsonBrowningTest

class VariableThresholdCMCTest: public ModelFitter{
    // write result header
    void writeHeader(FILE* fp) {
    };
    // fitting model
    int fit(VCFData& data) {
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
    };
}; // VariableThresholdCMCTest

class VariableThresholdFreqTest: public ModelFitter{
    // write result header
    void writeHeader(FILE* fp) {
    };
    // fitting model
    int fit(VCFData& data) {
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
    };
}; // VariableThresholdFreqTest

class SkatTest: public ModelFitter{
    // write result header
    void writeHeader(FILE* fp) {
    };
    // fitting model
    int fit(VCFData& data) {
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
    };
}; // SkatTest

#if 0

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

#endif


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
