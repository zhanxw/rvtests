#ifndef _VCFDATA_H_
#define _VCFDATA_H_

#include "OrderedMap.h"
#include "MathMatrix.h"
#include "VCFUtil.h"
#include "PlinkInputFile.h"
#include "PlinkOutputFile.h"

const static int FREQ_ALL = 1;
const static int FREQ_CONTORL_ONLY = 2;

/**
 * Hold genotype, phenotype, covariate data
 *      markerName (id), markerChrom, markerPos, markerRef, markerAlt, markerFreq
 *      peopleName
 *      marker2Idx, people2Idx
 * Read from VCF file
 * Read from PLINK format
 * Read external phenotype, covariate file
 */
class VCFData{
public:
    VCFData(){
        this->genotype = new Matrix;
        this->collapsedGenotype = new Matrix;
        this->phenotype = new Matrix;
        this->covariate = new Matrix;        
        this->phenotypeVec = new Vector;
        assert(this->genotype && this->phenotype && this->covariate);
    }
    ~VCFData(){
        if (this->genotype ) delete this->genotype ;
        if (this->collapsedGenotype) delete this->collapsedGenotype;
        if (this->phenotype) delete this->phenotype;
        if (this->covariate) delete this->covariate;
        if (this->phenotypeVec) delete this->phenotypeVec;

        this->genotype = NULL;
        this->collapsedGenotype = NULL;
        this->phenotype = NULL;
        this->covariate = NULL;
        this->phenotypeVec = NULL;
    }

    /**
     * will adjust people2idx and phenotype and covariate order
     * @param p: new order
     * @param droppedId: dropped people's id
     */
    void adjustPeopleOrder(std::vector< std::string>& p, std::vector<std::string>* droppedId){
        std::map<int, std::string> order; // will adjust internal order of phenotype and covariate.
        std::map<int, std::string>::iterator order_it;

        int nDropped = 0;
        for (int i = 0; i < p.size(); i++){
            if (!this->people2Idx.find(p[i])) { // VCF file has p[i], but not in existing phenotype or covariate.
                nDropped ++;
                if (droppedId)
                    droppedId->push_back(p[i]);
            } else{
                order[ order.size() ] = p[i];
            }
        }

        // adjust phenotype and covaraite order so they are the same as VCF file
        for (int i = 0; i < order.size() ; i++){
            //switch row i so i have people order[i]
            int j = this->people2Idx[order[i]];
            if (i == j) continue;
            else{
                if (this->phenotype->rows > j )
                    this->phenotype->SwapRows(i,j);
                if (this->covariate->rows > j )
                    this->covariate->SwapRows(i,j);
            }
        }
        fprintf(stdout, "%d samples from VCF files will not be included in study.\n", nDropped);
        this->people2Idx.clear();
        for (int i = 0; i < order.size(); i ++)
            this->people2Idx[order[i]] = i;
    };

    // adjust covariate and pheonotype 
    void addVCFHeader(VCFHeader* h, std::vector<std::string>* droppedId){
        assert(h);

        std::vector<std::string> v;
        h->getPeopleName(&v); // include all people
        if (this->people2Idx.size() == 0) {
            // don't have phenotype or covariate,
            // use VCF order to set people2idx file
            for (int i = 0; i < v.size(); i++)
                this->people2Idx[v[i]] = i;
            return;
        }

        // already have phenotype and covariate data
        assert(droppedId);
        this->adjustPeopleOrder(v, droppedId);
    };
        
    /**
     * NOTE: call addVCFHeader()  first
     */
    void addVCFRecord(VCFRecord& r){
        // make sure the phenotype, covariate and matched by peopleID
        // and match genotype peopleID to pheontoype and covaraite

        // add people genotype
        // 1. find markerName/row index
        VCFPeople& people = r.getPeople();
        std::string m = r.getID();
        if (m == ".") {
            m = r.getChrom();
            m += ':';
            m += toString(r.getPos());
        }
        int rowNum = -1;
        if (this->marker2Idx.find(m)) {
            rowNum = this->marker2Idx[m];
        } else{
            rowNum = this->genotype->rows;
            this->marker2Idx[m] = rowNum;
            this->genotype->Dimension(rowNum + 1, people.size());
        };
        

        VCFIndividual* indv;
        for (int i = 0; i < people.size(); i++) {
            indv = people[i];
            // find peopleName/col index.
            if (!this->people2Idx.find(indv->getName())) {
                continue;
            }
            int colNum = this->people2Idx[indv->getName()];

            // get GT index. if you are sure the index will not change, call this function only once!
            int GTidx = r.getFormatIndex("GT");
            if (GTidx >= 0) 
                (*this->genotype)[rowNum][colNum] = (*indv)[GTidx].getGenotype();
            else 
                (*this->genotype)[rowNum][colNum] = MISSING_GENOTYPE;
        }

    };
    // load data from plink format
    void loadPlink(const char* prefix){
        std::string p = prefix;
        if (p.size() == 0){
            fprintf(stderr, "Cannot open PLINK file %s.\n", prefix);
            return;
        }
        // load marker (BIM)
        this->loadMarkerFromBim( (p + ".bim").c_str());

        // load people (FAM)
        this->loadPeopleFromFam( (p + ".fam").c_str());
        
        // load bed (BED) into memory (??may use mmap() to shrink memory usage)
        // check magic word and snp major
        // read all rests into memory
        char magic[2];
        char mode; 
        FILE* fBed = fopen( (p + ".bed").c_str(), "rb");
        fread(magic, sizeof(char), 2, fBed);
        fread(&mode, sizeof(char), 1, fBed);
        if (magic[0] != 0x6c || magic[1] != 0x1b) {
            fprintf(stderr, "Cannot open BED file %s, corrupt magic word.\n", prefix);
            return;
        }

        // we reverse the two bits as defined in PLINK format, 
        // so we can process 2-bit at a time.
        const static unsigned char HOM_REF = 0x0;     //0b00;
        const static unsigned char HET = 0x2;         //0b10;
        const static unsigned char HOM_ALT = 0x3;     //0b11;
        const static unsigned char MISSING = 0x1;     //0b01;

        if (mode == 0x01) {
            // snp major mode
            unsigned char mask[] = { 0x3, 0xc, 0x30, 0xc0 }; //0b11, 0b1100, 0b110000, 0b11000000
            unsigned char c;
            int numPeople = this->people2Idx.size();
            int numMarker = this->marker2Idx.size();
            (*this->genotype).Dimension( numMarker, numPeople);
            for (int m = 0; m < numMarker; m++){
                for (int p = 0; p < numPeople; p++) {
                    int offset = p & (4 - 1);
                    if (offset == 0) {
                        fread(&c, sizeof(unsigned char), 1, fBed);
                    }
                    unsigned char geno = (c & mask[offset]) >> (offset << 1);
                    switch (geno){
                    case HOM_REF:
                        (*this->genotype)[m][p] = 0;
                        break;
                    case HET:
                        (*this->genotype)[m][p] = 1;
                        break;
                    case HOM_ALT:
                        (*this->genotype)[m][p] = 2;
                        break;
                    case MISSING:
                        (*this->genotype)[m][p] = MISSING_GENOTYPE;
                        break;
                    default:
                        REPORT("Read PLINK genotype error!\n");
                        break;
                    };
                }
            }
        } else if (mode == 0x00) {
            // people major mode
            // TODO
            fprintf(stderr, "Cannot open BED file %s, work to be DONE.\n", prefix);
        } else {
            fprintf(stderr, "Cannot open BED file %s, unrecognized major mode.\n", prefix);
        };
        fclose(fBed);
    };

    // PLINK requires phenotype on column 3.
    int loadPlinkPhenotype(const char* fn) {
        return this->loadPhenotypeByCol(fn, 3);
    };
    // col: 1-based column 
    // return: num of people whose phenotype that are not set
    // return 0: success ; -n: n individual with missing phenotype
    int loadPhenotypeByCol(const char* fn, int col) { // by default PLINK use 3rd column as phenotype
        // if this->genotype already have people2Idx, then make phenotype people label match (also in order)
        // otherwise, just set people2Idx like the phenotype.
        if (this->people2Idx.size()) {
            OrderedMap<std::string, int> pName;
            int ret = this->readPlinkPhenotypeSkipMissing(fn, 0,
                                                          this->phenotype,
                                                          &pName,
                                                          &this->phenotype2Idx);
            if (ret) {
                fprintf(stderr, "%d phenotype are skipped as their pheontype are not numeric.\n", ret);
            } 
            int skip1, skip2;
            this->matchData(this->genotype, this->people2Idx, &skip1, 
                            this->phenotype, pName, &skip2);
            if (skip1) {
                fprintf(stderr, "%d individuals' genotypes are dropped because of mismatch.\n", skip1);
            } 
            if (skip2) {
                fprintf(stderr, "%d individuals' phenotype are dropped because of mismatch.\n", skip2);
            }
            return ret;
        } else {

            int ret = this->readPlinkPhenotypeSkipMissing(fn, 
                                                          0,
                                                          this->phenotype,
                                                          &this->people2Idx,
                                                          &this->phenotype2Idx);
            if (ret){
                fprintf(stderr, "%d phenotype are skipped as their pheontype are not numeric.\n", ret);
            }
            return ret;
        }
    };
    int readPlinkPhenotypeSkipMissing(const char* fn, const char* selectedCol, 
                                      Matrix* data, 
                                      OrderedMap<std::string, int> * rowName,
                                      OrderedMap<std::string, int> * colName);

    // abover @param threshold will set to 1.0, or 0.0
    void dichotomizedPhenotype(double threshold) {
        for (int i = 0; i < phenotypeVec->Length(); i++){
            if ((*this->phenotype)[i][0] > threshold) {
                (*this->phenotype)[i][0] = 1.0;
                (*phenotypeVec)[i] = 1.0;
            } else {
                (*this->phenotype)[i][0] = 0.0;
                (*phenotypeVec)[i] = 0.0;
            }
        }
    };

    Vector* extractPhenotype(){
        if (this->phenotypeVec->Length()) return this->phenotypeVec;

        this->phenotypeVec->Dimension(this->people2Idx.size());
        for (int i = 0; i < phenotypeVec->Length(); i++){
            (*phenotypeVec)[i] = (*this->phenotype)[i][0];
        }
        return this->phenotypeVec;
    };
    // we only load covariate for people that appeared in genotype
    // please use numeric covariate only
    // return: num of people whose covariate that are not set
    // return 0: success
    int loadCovariate(const char* fn){
        if (!fn || strlen(fn) == 0){
            return 0;
        }
        double defaultMissingCovariate = 0.0;
        // if this->genotype already have people2Idx, then make covariate people label match (also in order)
        // otherwise, just set people2Idx like the covariate.
        if (this->people2Idx.size()) {
            OrderedMap<std::string, int> pName;
            int ret = this->readPlinkTable(fn, 
                                           this->covariate,
                                           &pName, 
                                           &this->covariate2Idx,
                                           defaultMissingCovariate);
            int skip1, skip2;
            this->matchData(this->genotype, this->people2Idx, &skip1,
                            this->covariate, pName, &skip2);
            if (skip1) {
                fprintf(stderr, "%d individuals' genotypes are dropped because of mismatch.\n", skip1);
            } 
            if (skip2) {
                fprintf(stderr, "%d individuals' covariates are dropped because of mismatch.\n", skip2);
            }
            return ret;
        } else {
            int ret = this->readPlinkTable(fn, 
                                           this->covariate,
                                           &this->people2Idx, 
                                           &this->covariate2Idx,
                                           defaultMissingCovariate);
            if (ret){
                fprintf(stderr, "%d covariates are set to missing value %.2f.\n", ret, defaultMissingCovariate);
            }
            return ret;
        }
    };
    // loop through rowLabel1, and switch the command rowLabel to the toppest position
    void matchData(Matrix* m1, OrderedMap<std::string, int>& rowLabel1, int* skip1,
                   Matrix* m2, OrderedMap<std::string, int>& rowLabel2, int* skip2) {
        assert (skip1 && skip2);

        OrderedMap<std::string, int> newLabel;
        int nextRow1 = 0;
        int nextRow2 = 0;

        for (unsigned int i = 0; i < rowLabel1.size() ; i++){
            const std::string& key1 = rowLabel1.keyAt(i);
            if (rowLabel2.find(key1)){
                // both matrix have the row
                int j = rowLabel2[key1];
                newLabel[key1] = nextRow1;
                if (nextRow1 != i)
                    m1->SwapRows(nextRow1, i);
                nextRow1++;
                if (nextRow2 != i) 
                    m2->SwapRows(nextRow2, j);
                nextRow2++;
            } else {
                // only m1 have that row
                continue;
            }
        }
        *skip1 = rowLabel1.size() - nextRow1;
        *skip2 = rowLabel2.size() - nextRow2;
        rowLabel1 = newLabel;
        rowLabel2 = newLabel;
    };

    /**
     * read @param fn into @param data from R-readable format
     * REQURIE:
     *   @param rowName and @param colName are both treated like string
     *   @param data should be integer/float number
     *   @param defaultValue, when data cannot be converted to double, use this value instead.
     */
    int readTable(const char* fn,
                  Matrix* data, 
                  OrderedMap<std::string, int> * rowName,
                  OrderedMap<std::string, int> * colName,
                  std::string* upperLeftName,
                  double defaultValue);
    /**
     * similar to readTable, but skip the first column and
     * require header file have the same number of fields as other lines.
     */
    int readPlinkTable(const char* fn,
                       Matrix* data, 
                       OrderedMap<std::string, int> * rowName,
                       OrderedMap<std::string, int> * colName,
                       double defaultValue);
    // read genotypes.
    // required format:
    //   have header: PeopleID MarkerName[0] MarkerName[1]...
    //   don't have row names
    //   genotype are people x marker
    void readGenotypeFromR(const char* fn){
        this->people2Idx.clear();
        this->marker2Idx.clear();

        if (!fn || strlen(fn) == 0) {
            fprintf(stderr, "Cannot read genotypes from %s, as it is empty.\n", fn);
            return;
        }
        
        LineReader lr(fn);
        std::vector <std::string> fd;
        int lineNo = 0;
        while (lr.readLineBySep(&fd, "\t ")){
            if (!lineNo) {// header line
                for (int i = 1; i < fd.size(); i++) {
                    if (this->marker2Idx.find(fd[i])) {
                        fprintf(stderr, "Duplicate marker %s.\n", fd[i].c_str());
                    }
                    this->marker2Idx[fd[i]] = (i-1);
                };
                fprintf(stdout, "%d marker(s) loaded.\n", this->marker2Idx.size());
            } else { // body lines
                this->people2Idx[fd[0]] = (lineNo - 1);
                this->genotype->Dimension( this->marker2Idx.size(), this->people2Idx.size());
                for (unsigned int  m = 1; m < fd.size(); m++) {
                    (*this->genotype) [ (lineNo - 1)] [ int(m)] = atoi(fd[m]);
                }
            }
            lineNo ++;
        };
    };

// outputs:
// prefix + ".geno": raw genotype
// prefix + ".cgeno": raw genotype
// prefix + ".cov": raw genotype
// prefix + ".pheno": raw genotype
    void writeRawData(const char* prefix);
// write genotype in R format:
//   have header: PeopleID MarkerName[0] MarkerName[1]...
//   don't have row names, e.g. 1, 2, 3,...
//   genotype are people x marker
    void writeGenotype(const char* fn);
    // void writeCollapsedGenotype( const char* fn);
    void writeCovariate(const char* fn);
    void writePhenotype(const char* fn);
    void writeTable(const char* fn,
                    Matrix* data, OrderedMap<std::string, int> & rowName,
                    OrderedMap<std::string, int> & colName,
                    const char* upperLeftName);

    
    // type: FREQ_ALL, FREQ_CONTORL_ONLY
    // will also update:
    //    markerMAC and markerTotoalAllele ( for above type)
    void calculateFrequency(const int type){
        if (type != FREQ_ALL && type != FREQ_CONTORL_ONLY) {
            fprintf(stderr, "Unsupported frequency calculation type!\n");
            abort();
        };
        int numMarker = this->genotype->rows;
        int numPeople = this->genotype->cols;
        this->markerFreq.resize(numMarker);
        this->markerMAC.resize(numMarker);
        this->markerTotalAllele.resize(numMarker);

        for (int m = 0; m < numMarker; m++){
            for (int p = 0; p < numPeople; p++) {
                if ( ! (
                         (type == FREQ_ALL) ||
                         (type == FREQ_CONTORL_ONLY && ((*this->phenotype)[p][0] - 1.0 < 1e-10))
                         )
                    )
                    continue;

                int geno = (*this->genotype)[m][p];
                if (geno >= 0){
                    this->markerTotalAllele[m] ++;
                    if (geno != 0){
                        this->markerMAC[m] += geno;
                    }
                }
            }
            if (this->markerTotalAllele[m] > 0)
                this->markerFreq[m] = (double)(this->markerMAC[m]) / this->markerTotalAllele[m];
            else 
                this->markerFreq[m] = 0.0;
        };
    };

    bool isCaseControlPhenotype(){
        Vector* v = this->extractPhenotype();
        // case control are coded as 0, 1, or missing;
        // if there are non-integer phenotype, then quantative trait
        // if case-control but not code as 0, 1, or -9, give warning!
        
        int doubleCnt = 0;
        int wrongCoding = 0;
        for (int i = 0; i < v->Length(); i++){
            double res = (*v)[i] - round( (*v)[i]);
            if (res > 1e-3 || res < -1e-3) doubleCnt ++;
            int code = (int) (*v)[i];
            if (code != 0 || code != 1 || code != -9) {
                wrongCoding ++;
            }
        }
        if (doubleCnt > 0) {return false;};
        if (wrongCoding > 0) {
            fprintf(stderr, "Case control coding should be 0, 1, or -9, but wrong coding detected.\n");
        }
        return true;
    };

    void dumpSize() {
        printf("genotype: %d x %d ;\n", genotype->rows, genotype->cols);
        printf("collapsedGenotype: %d x %d ;\n", collapsedGenotype->rows, collapsedGenotype->cols);
        printf("phenotype: %d x %d ;\n", phenotype->rows, phenotype->cols);
        printf("covariate: %d x %d ;\n", covariate->rows, covariate->cols);
        printf("people2Idx: %d ;\n", (int)people2Idx.size());
        printf("marker2Idx: %d ;\n", (int)marker2Idx.size());
    };

    /* // use friend class to save codes... */
    /* friend class Collapsor; */
    /* Matrix* getGeno() {return this->genotype;}; */
    /* Matrix* getPheno() {return this->phenotype;}; */
    /* Matrix* getCov() {return this->covariate;}; */
    /* OrderedMap<std::string, int>* getPeople2Idx() {return &this->people2Idx;}; */
    /* OrderedMap<std::string, int>* getMarker2Idx() {return &this->marker2Idx;}; */
private:
    void loadMarkerFromBim(const char* fn){
        std::vector<std::string> fd;
        LineReader lr(fn);
        int lineNo = 0;
        while(lr.readLineBySep(&fd, "\t ")){
            if (fd.size() != 6){
                fprintf(stderr, "BIM file %s corrupted.\n", fn);
                return;
            }
            this->marker2Idx[fd[1]] = lineNo++;
        };
        
        if (this->marker2Idx.size() != lineNo) {
            fprintf(stderr, "Duplicate markers in %s.\n", fn );
        }
    }
    void loadPeopleFromFam(const char* fn){
        std::vector<std::string> fd;
        LineReader f(fn);
        int lineNo = 0;
        Vector pheno;
        while(f.readLineBySep(&fd, "\t ")){
            if (fd.size() != 6){
                fprintf(stderr, "FAM file %s corrupted.\n", fn);
                return;
            }
            this->people2Idx[fd[1]] = lineNo++;
            pheno.Push(atoi(fd[5]));
        };
        if (this->people2Idx.size() != lineNo){
            fprintf(stderr, "Dupliate people name in %s.\n", fn);
        };
        // this->numPeople = this->people2Idx.size();
        (*this->phenotype).Dimension(this->people2Idx.size(), 1);
        
        for (int i = 0; i < pheno.Length(); i++){
            (*this->phenotype)[i][0] = pheno[i];
        }
    }
  public:
// member variables here are made public
// so that access is easier.
  public:
    Matrix* genotype;     // marker x people
    Matrix* collapsedGenotype; // people x marker
    Matrix* covariate;    // people x cov
    Matrix* phenotype;    // people x phenotypes
    Vector* phenotypeVec; // pheotype to be use in various statistical models. 
                          // here use it to speed up data loading.

    /* std::vector<std::string> markerName; */
    /* std::vector<int> markerChrom; */
    /* std::vector<int> markerPos; */
    /* std::vector<std::string> markerRef; */
    /* std::vector<std::string> markerAlt; */
    std::vector<double> markerFreq;
    std::vector<int> markerMAC;          // minor(alternative) allele count 
    std::vector<int> markerTotalAllele;  // Total number of allele

    std::vector<double> collapsedMarkerFreq;
    std::vector<int>    collapsedMarkerMAC;          // minor(alternative) allele count 
    std::vector<int>    collapsedMarkerTotalAllele;  // Total number of allele


    OrderedMap<std::string, int> people2Idx; // peopleID -> idx in this->genotype
    // int numPeople; 

    OrderedMap<std::string, int> marker2Idx; // markerName -> idx in this->genotype
    // int numMarker;

    OrderedMap<std::string, int> covariate2Idx; // covariate name -> idx in this->covariate
    OrderedMap<std::string, int> phenotype2Idx; // phenotype -> idx in this->collapsedGenotype
    // OrderedMap<std::string, int> set2Idx; // collapsedSetName -> idx in this->collapsedGenotype

}; // end VCFData

#endif /* _VCFDATA_H_ */
