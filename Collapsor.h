#ifndef _COLLAPSOR_H_
#define _COLLAPSOR_H_

#include "VCFData.h"

class Collapsor{
  public:
    void setSetFileName(const char* fn){
        this->setContent.clear();
        this->setName.clear();

        if (!fn || strlen(fn))
            return;

        LineReader lr(fn);
        std::vector<std::string> fd;
        int lineNo = 0;
        while(lr.readLineBySep(&fd, " \t")){
            lineNo ++;
            if (fd.size() < 3){
                fprintf(stderr, "Cannot recognized format on line %d. \n", lineNo);
                continue;
            }
            if (fd[1] == "RANGE"){
                RangeList rl;
                for (int i = 2; i < fd.size(); i++){
                    rl.addRangeList(fd[i].c_str());
                }
                this->setContent.push_back (rl);
                this->setName.push_back(fd[0]);

            } else if (fd[1] == "MARKER"){
#pragma message "TODO"
                fprintf(stderr, "Unsupported for now \n", lineNo);
            } else {
                fprintf(stderr, "Cannot regconized keyword on line %d. \n", lineNo);
            }
        };
        this->setIndex = -1;  //reset setIndex
    };
    void setCollapsingStrategy(const int strategy){
#pragma message "Check if compatible, e.g NAIVE and set file are not compatible"
        this->collapsingStrategy = strategy;
    };
    /** iterate  all sets
     *  @return false: when all sets are went over
     */
    bool iterateSet(VCFInputFile& vin, VCFData* data){
        setIndex ++;
        // TODO: match individuals
#pragma messge "Handle sample matching program"
        //data->addVCFHeader(vin.getVCFHeader());
        data->addVCFHeader(vin.getVCFHeader());

        if (this->setName.size() == 0) {
            // iterate every marker
            if (vin.readRecord()){
                VCFRecord& record = vin.getVCFRecord();
                data->addVCFRecord(record);

                this->currentSetName =record.getID();
                if (this->currentSetName == ".") {
                    this->currentSetName = record.getChrom();
                    this->currentSetName += toString(record.getPos());
                }
                return true;
            } else{
                return false;
            }
        }

        // check boundary condition
        if (setIndex == this->setName.size()) return false;

        // load one set
        this->currentSetName = this->setName[this->setIndex];
        RangeList& rl = this->setContent[this->setIndex];

        // add genotypes within set to this-> genotype
        for (int i = 0; i < rl.size(); i++ ){
            vin.setRangeList(rl);
            while(vin.readRecord()){
                VCFRecord& record = vin.getVCFRecord();
                data->addVCFRecord(record);
            }
        };
    };

    std::string& getCurrentSetName() {
        return this->currentSetName;
    };

    /**
     * Collapse data->genotype to this->collapsedGenotype
     * @param param can be used to specify additional parameters.
     */
    bool collapseGenotype(VCFData* data, void* param){
        // note: missing data are handled different
        // see each collapsing method implementation for details.
        switch (this->collapsingStrategy) {
        case NAIVE:
            this->naiveCollapse(data);
            break;
        case CMC:
            this->cmcCollapse(data);
            break;
        case ZEGGINI:
            this->zegginiCollapse(data);
            break;
        case MADSON_BROWNING:
            this->madsonbrowningCollapse(data);
            break;
        case PROGRESSIVE:
            this->progressiveCollapse(data, param);
            break;
        case UNDEFINED:
            fprintf(stderr, "Please call setCollapsingStrategy().\n");
            abort();
        default:
            fprintf(stderr, "Unrecognized collapsing method!\n");
            return false;
        };
        // sanity check collapsing results
        assert(data->collapsedGenotype->cols == data->set2Idx.size());
    };
    bool outputHeader(FILE* fp){
        fprintf(fp, "%s", this->currentSetName.c_str());
    };
  public:
    void naiveCollapse(VCFData* d){
        int numPeople = d->people2Idx.size();
        int numMarker = d->marker2Idx.size();
        d->collapsedGenotype->Dimension(numPeople, numMarker);
        for (int p = 0; p < numPeople; p++)
            for (int m = 0; m < numMarker; m++)
                (*d->collapsedGenotype)[p][m] = (*d->genotype)[m][p];
        d->set2Idx = d->marker2Idx;
    };
    void cmcCollapse(VCFData* d){
        int numPeople = d->people2Idx.size();
        int numMarker = d->marker2Idx.size();
        d->collapsedGenotype->Dimension(numPeople, 1);
        d->collapsedGenotype->Zero();
        for (int p = 0; p < numPeople; p++){
            for (int m = 0; m < numMarker; m++) {
                int g = (int)((*d->genotype)[m][p]);
                if (g > 0) {
                    (*d->collapsedGenotype)[p][0] = 1.0;
                    break;
                }
            };
        };
    };
    void zegginiCollapse(VCFData* d){
        int numPeople = d->people2Idx.size();
        int numMarker = d->marker2Idx.size();
        d->collapsedGenotype->Dimension(numPeople, 1);
        d->collapsedGenotype->Zero();
        for (int p = 0; p < numPeople; p++){
            for (int m = 0; m < numMarker; m++) {
                int g = (*d->genotype)[m][p];
                if (g > 0) { // genotype is non-reference
                    (*d->collapsedGenotype)[p][0] += 1.0;
                    break;
                }
            };
        };
    };
    void madsonbrowningCollapse(VCFData* d){
        d->calculateFrequency(VCFData::FREQ_CONTORL_ONLY);

        int numPeople = d->people2Idx.size();
        int numMarker = d->marker2Idx.size();
        d->collapsedGenotype->Dimension(numPeople, 1);
        d->collapsedGenotype->Zero();

        for (int m = 0; m < numMarker; m++) {
            // up weight by control freuqncey  1/p/(1-p)
            double weight = 0.0;
            if (d->markerFreq[m] >= 1e-6){
                weight = 1.0 / d->markerFreq[m] / (1.0 - d->markerFreq[m]);
            } else{
                continue;
            }

            // calculate burden score at marker m
            for (int p = 0; p < numPeople; p++){
                double g = (*d->genotype)[m][p];
                if (g > 0) {
                    (*d->collapsedGenotype)[p][0] += g * weight;
                    break;
                }
            };
        };
    };
    void progressiveCollapse(VCFData* d, void* param){
        int col = *(int*)param;
        int collapsingStrategy = * ( (int*)(param) + 1);
        if (collapsingStrategy != CMC && collapsingStrategy != MADSON_BROWNING){
            fprintf(stderr, "Unknow progressive collapsing method: %d.\n", collapsingStrategy);
            abort();
        }

        int numPeople = d->people2Idx.size();
        int numMarker = d->marker2Idx.size();
        if (col < 0){
            d->calculateFrequency(VCFData::FREQ_CONTORL_ONLY);
            d->collapsedGenotype->Dimension(numPeople, 1);
            d->collapsedGenotype->Zero();
            if (collapsingStrategy = MADSON_BROWNING){
                d->calculateFrequency(VCFData::FREQ_CONTORL_ONLY);
            }
            return;
        }
        switch(collapsingStrategy){
        case CMC:
            for (int p = 0; p < numPeople; p++) {
                int g = (int)((*d->genotype)[col][p]);
                if (g > 0) {
                    (*d->collapsedGenotype)[p][0] = 1.0;
                    break;
                }
            };
            break;
        case MADSON_BROWNING:
        {
            double weight = 0.0;
            if (d->markerFreq[col] >= 1e-6){
                weight = 1.0 / d->markerFreq[col] / (1.0 - d->markerFreq[col]);
            } else{
                return;
            }
            // calculate burden score at marker m
            for (int p = 0; p < numPeople; p++){
                double g = (*d->genotype)[col][p];
                if (g > 0) {
                    (*d->collapsedGenotype)[p][0] += g * weight;
                    break;
                }
            };
        }
        break;
        default:
            fprintf(stderr, "Wrong progressive collapsing method: %d.\n", collapsingStrategy);
            abort();
            break;
        }
    };

  public:
    static const int UNDEFINED = -1;
    static const int NAIVE = 0;
    static const int CMC = 1;
    static const int ZEGGINI = 2;
    static const int MADSON_BROWNING = 3;
    static const int PROGRESSIVE = 4;

  private:
    Matrix stagedGenotype; // marker x people
    int collapsingStrategy;
    OrderedMap<std::string, int>* people2Idx;
    OrderedMap<std::string, int>* marker2Idx;

    std::string currentSetName;

    int setIndex;                       // record current visited index of set
    std::vector< std::string> setName;
    std::vector<RangeList> setContent;
};

#endif /* _COLLAPSOR_H_ */
