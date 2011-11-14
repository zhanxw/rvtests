#ifndef _VCFDATA_H_
#define _VCFDATA_H_

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
        this->phenotype = new Matrix;
        this->covariate = new Matrix;        
        assert(this->genotype && this->phenotype && this->covariate);
    }
    ~VCFData(){
        if (this->genotype ) delete this->genotype ;
        if (this->phenotype) delete this->phenotype;
        if (this->covariate) delete this->covariate;
    }
    void loadVCFHeader(VCFHeader* h){
        std::vector<std::string> p;
        h->getPeopleName(&p);
        this->numPeople = p.size();
        for (int i = 0; i < p.size(); i++){
            this->people2Idx[p[i]] = i;
        };
    };
    void loadVCFRecord(VCFRecord* r){
        //add one marker for all people
        FATAL("Not ready at this moment. Will provide soon.\n");
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
            if (this->people2Idx.count(pname) == 0) {
                fprintf(stderr, "%s does not exist yet.\n", pname.c_str());
                continue;
            }
            int idx = this->people2Idx[pname];
            processed.insert(pname);
            (*this->phenotype)[idx][0] = atof(fd[col].c_str());
        };
        return (this->people2Idx.size() - processed.size());
    };
    // we only load covariate for people that appeared in genotype
    // please use numeric covariate only
    // return: num of people whose covariate that are not set
    // return 0: success
    int loadCovariate(const char* fn){
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
            if (this->people2Idx.count(pname) == 0) {
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
    void writeGenotypeToR(const char* fn){
        FileWriter fw(fn);
        // header
        for (int i = 0; i < markerName.size(); i++){
            if (i) fw.write(" ");
            if (markerName[i].size() == 0){
                fw.write("\".\"");
            }else {
                fw.write("\"");
                fw.write(markerName[i].c_str());
                fw.write("\"");
            }
        };
        if (markerName.size()) fw.write("\n");

        // content
        for (int m = 0; m < this->numMarker; m++) {
            fw.printf("%d", m); // line no
            for (int p = 0; p < this->numPeople; p++ ){
                fw.printf(" %d", (int)((*this->genotype)[m][p]));
            }
            fw.write("\n");
        }
    };

    // use friend class to save codes...
    friend class Collapsor;
    Matrix* getGeno() {return this->genotype;};
    Matrix* getPheno() {return this->phenotype;};
    Matrix* getCov() {return this->covariate;};
    std::map<std::string, int>* getPeople2Idx() {return &this->people2Idx;};
    std::map<std::string, int>* getMarker2Idx() {return &this->marker2Idx;};
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
        this->numMarker = marker2Idx.size();
        
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
        this->numPeople = this->people2Idx.size();
        (*this->phenotype).Dimension(numPeople, 1);
        
        for (int i = 0; i < pheno.Length(); i++){
            (*this->phenotype)[i][0] = pheno[i];
        }
    }

private:
    Matrix* genotype; // marker x people
    Matrix* covariate; // people x cov
    Matrix* phenotype; // people x phenotypes
    
    std::vector<std::string> markerName;
    std::vector<int> markerChrom;
    std::vector<int> markerPos;
    std::vector<std::string> markerRef;
    std::vector<std::string> markerAlt;
    std::vector<double> markerFreq;
    std::vector<int> markerCount;

    std::map<std::string, int> people2Idx;
    int numPeople;

    std::map<std::string, int> marker2Idx;
    int numMarker;
}; // end VCFData

#endif /* _VCFDATA_H_ */
