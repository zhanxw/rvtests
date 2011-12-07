#ifndef _VCFDATA_H_
#define _VCFDATA_H_

#include "OrderedMap.h"
#include "MathMatrix.h"

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
        assert(this->genotype && this->phenotype && this->covariate);
    }
    ~VCFData(){
        if (this->genotype ) delete this->genotype ;
        if (this->collapsedGenotype) delete this->collapsedGenotype;
        if (this->phenotype) delete this->phenotype;
        if (this->covariate) delete this->covariate;
        this->genotype = NULL;
        this->collapsedGenotype = NULL;
        this->phenotype = NULL;
        this->covariate = NULL;
    }
    // adjust covariate and pheonotype
    // so the data will be matched later on.
    void addVCFHeader(VCFHeader* h){
        std::vector<std::string> p;
        h->getPeopleName(&p);
        for (int i = 0; i < p.size(); i++){
            this->people2Idx[p[i]] = i;
        };
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
    // col: 1-based column 
    // return: num of people whose phenotype that are not set
    // return 0: success
    int overwritePhenotype(const char* fn, int col) { // by default PLINK use 3rd column as phenotype
        if (col <= 0) {
            fprintf(stderr, "col should be larger than 0.\n");
            return -1;
        }

        col --; // get 0-based column
        LineReader lr(fn);
        std::vector<std::string> fd;
        std::set<std::string> processed;
        while(lr.readLineBySep(&fd, "\t ")){
            if (fd.size() < 3 || col >= fd.size() ) {
                fprintf(stderr, "Insufficient columns in %s.\n", fn);
                continue;
            }
            const std::string& pname = fd[1];
            if (!this->people2Idx.find(pname) ) {
                fprintf(stderr, "%s does not exist yet.\n", pname.c_str());
                continue;
            }
            int idx = this->people2Idx[pname];
            processed.insert(pname);
            (*this->phenotype)[idx][0] = atof(fd[col].c_str());
        };
        return (this->people2Idx.size() - processed.size());
    };
    void extractPhenotype(Vector* v){
        v->Dimension(this->people2Idx.size());
        for (int i = 0; i < v->Length(); i++){
            (*v)[i] = (*this->phenotype)[i][0];
        }
    };
    // we only load covariate for people that appeared in genotype
    // please use numeric covariate only
    // return: num of people whose covariate that are not set
    // return 0: success
    int readCovariate(const char* fn){
        // we cannot use phenotype->Read(FILE*) ,
        // since in that case header line is required!
        LineReader lr(fn);
        std::vector<std::string> fd;
        std::set<std::string> processed;
        int fdLen = -1;
        while(lr.readLineBySep(&fd, "\t ")){
            if (fd.size() < 3 ) {
                fprintf(stderr, "Insufficient columns in %s.\n", fn);
                continue;
            }
            if (fdLen > 0 && fdLen != fd.size()) {
                fprintf(stderr, "Inconsistent columns in %s.\n", fn);
                continue;
            };
            if (fdLen < 0) {// first line
                fdLen = fd.size(); 
                // two leading columns are Family ID (unsed) and Individual ID (used).
                this->covariate->Dimension(this->people2Idx.size(), fd.size() - 2);
            }
            const std::string& pname = fd[1];
            if (!this->people2Idx.find(pname)) {
                fprintf(stderr, "%s does not exist yet.\n", pname.c_str());
                continue;
            }
            int idx = this->people2Idx[pname];
            processed.insert(pname);
            for (int col = 2; col < fd.size(); col ++) 
                (*this->covariate)[idx][col - 2] = atof(fd[col].c_str());
        };
        return (this->people2Idx.size() - processed.size());
    };
    int readTabularFile(const char* fn, Matrix* m, OrderedMap<std::string, int>* rowLabel, OrderedMap<std::string, int>* colLabel){
        fprintf(stderr, "TODO\n");
        return 0;
    };
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
    void writeRawData(const char* prefix){
        std::string p = prefix;
        if (p == "") {
            p = "rvtest.raw";
        }
        this->writeGenotype( (p + ".geno").c_str());
        this->writeCollapsedGenotype( (p + ".cgeno").c_str());
        this->writeCovariate( (p + ".cov").c_str());
        this->writePhenotype( (p + ".pheno").c_str());
    };

    // write genotype in R format:
    //   have header: PeopleID MarkerName[0] MarkerName[1]...
    //   don't have row names
    //   genotype are people x marker
    void writeGenotype(const char* fn){
        this->writeTable(fn, this->genotype, this->marker2Idx, this->people2Idx, "MarkerName");
    };
    void writeCollapsedGenotype( const char* fn){
        this->writeTable(fn, this->collapsedGenotype, this->people2Idx, this->set2Idx, "PeopleID");
    };
    /** write covariate to file, format is as following:
     *  header line: PeopleName CovName1 CovName2 ...
     *  content line: P1 1.0 2.0 ...
     */
    void writeCovariate(const char* fn) {
        this->writeTable(fn, this->covariate, this->people2Idx, this->covariate2Idx, "PeopleID");
    };
    void writePhenotype(const char* fn) {
        this->writeTable(fn, this->phenotype, this->people2Idx, this->phenotype2Idx, "PeopleID");
    };

    /**
     */
    void writeTable(const char* fn,
                    Matrix* data, OrderedMap<std::string, int> & rowName, 
                    OrderedMap<std::string, int> & colName, 
                    const char* upperLeftName) {
        if (!data || data->rows == 0 || data->cols == 0)
            return;

        if (rowName.size() != data->rows) {
            fprintf(stderr, "Row number does not match!");
            return;
        }
        if (colName.size() != data->cols) {
            fprintf(stderr, "Col number does not match!");
            return;
        }

        FileWriter fw(fn);
        // header
        fw.write(upperLeftName);
        for (int i = 0; i < colName.size(); i++){
            fw.write("\t");
            if (colName.size() == 0){
                fw.write("\".\"");
            }else {
                fw.write("\"");
                fw.write(colName.keyAt(i).c_str());
                fw.write("\"");
            }
        };
        fw.write("\n");

        // content
        int numCol = colName.size();
        int numRow = rowName.size();
        for (int r = 0; r < numRow; r++ ){
            fw.write(rowName.keyAt(r).c_str());
            for (int c = 0; c < numCol; c++) {
                fw.printf("\t%d", (int)((*data)[r][c]));
            }
            fw.write("\n");
        }
        fw.close();
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
    
    void calcFreqAndCount() {
        this->markerFreq.clear();
        this->markerCount.clear();
        this->markerTotalAllele.clear();
            
        for (int i = 0; i < genotype->rows; i++) {
            int count = 0;
            int totalAllele = 0;
            for (int j = 0; j < genotype->cols; j++) {
                if (genotype < 0) continue;
                count += (*genotype)[i][j];
                totalAllele += 2;
            }
            this->markerCount.push_back(count);
            this->markerTotalAllele.push_back(totalAllele);
            if (totalAllele == 0)
                this->markerFreq.push_back(-1);
            else
                this->markerFreq.push_back( (double)(count) / totalAllele);
        }
    };

// member variables here are made public
// so that access is easier.
  public:
    Matrix* genotype; // marker x people
    Matrix* collapsedGenotype;
    Matrix* covariate; // people x cov
    Matrix* phenotype; // people x phenotypes
    
    /* std::vector<std::string> markerName; */
    /* std::vector<int> markerChrom; */
    /* std::vector<int> markerPos; */
    /* std::vector<std::string> markerRef; */
    /* std::vector<std::string> markerAlt; */
    std::vector<double> markerFreq;
    std::vector<int> markerCount; // Alternative allele count
    std::vector<int> markerTotalAllele; // Total number of allele

    OrderedMap<std::string, int> people2Idx; // peopleID -> idx in this->genotype
    // int numPeople; 

    OrderedMap<std::string, int> marker2Idx; // markerName -> idx in this->genotype
    // int numMarker;

    OrderedMap<std::string, int> covariate2Idx; // covariate name -> idx in this->covariate
    OrderedMap<std::string, int> phenotype2Idx; // phenotype -> idx in this->collapsedGenotype
    OrderedMap<std::string, int> set2Idx; // collapsedSetName -> idx in this->collapsedGenotype

}; // end VCFData

#endif /* _VCFDATA_H_ */
