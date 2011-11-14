#ifndef _VCFUTIL_H_
#define _VCFUTIL_H_

#include "Exception.h"
#include "RangeList.h"
#include "Utils.h"

#include "libsrc/MathMatrix.h"

class VCFHeader{
  public:
    void push_back(const std::string& s){
        this->data.push_back(s);
    }
    void getPeopleName(std::vector<std::string>* p) const{
        if (!p) return;
        if (this->data.size() < 1) return;
        const std::string ln = this->data[this->data.size() - 1];
        
        std::vector<std::string> fd;
        stringTokenize(ln, "\t", &fd);
        if (fd.size() < 10) return;
        p->clear();
        for (int i = 9; i < fd.size() ; i++) {
            p->push_back(fd[i]);
        }
    };
    int size() const {return this->data.size();}
    std::string operator[] (int n) {return this->data[n];};
    const std::string operator[] (int n) const {return this->data[n];};
    std::string at(int n) { return this->data.at(n);};
    const std::string at(int n) const {return this->data.at(n);};
  private:
    std::vector<std::string> data;
};

#define MISSING_GENOTYPE -9
/**
 * the versatile format to store value for VCF file.
 */
class VCFValue{
public:
    int beg; // inclusive
    int end; // exclusive, and beg <= end
    const char* line;
public:
    int toInt() const{ 
        return atoi(line+beg);
    };
    void toInt(int* i) const {
        *i = atoi(line+beg);
    };
    double toDouble() const {
        return atoi(line+beg);        
    };
    void toDouble(double* d) const {
        *d = strtod(line+beg, 0);
    };
    /* std::string toStr() const {  */
    /*     std::string s; */
    /*     for (int i = beg; i < end; i++){ */
    /*         s.push_back(line[i]); */
    /*     } */
    /*     return s; */
    /* }; */
    const char* toStr() {
        this->retStr.clear();
        for (int i = beg; i < end; i++){
            this->retStr.push_back(line[i]);
        }
        return (this->retStr.c_str());
    };
    void toStr(std::string* s) const { 
        s->clear();
        for (int i = beg; i < end; i++){
            s->push_back(line[i]);
        }
    };
public:
    int getGenotype(){
        int g = 0;
        int p = beg;
        if (line[p] == '.') 
            return MISSING_GENOTYPE;
        if (line[p] < '0') 
            REPORT("Wrong genotype detected. [1]");
        else 
            g += line[p] - '0';
        
        p ++ ;
        if (p == end) 
            return g;
        
        p ++;
        if (p == end)
            REPORT("Wrong genotype length = 2");
        if (line[p] == '.')
            return MISSING_GENOTYPE;
        if (line[p] < '0')
            REPORT("Wrong genotype detected. [2]");
        else 
            g += line[p] - '0';
        
        return g;
    };
    int getAllele1(){
        int g = 0;
        int p = beg;
        if (line[p] == '.') 
            return MISSING_GENOTYPE;
        if (line[p] < '0') 
            REPORT("Wrong genotype detected. [1]");
        else 
            g += line[p] - '0';
        return g;
    };
    int getAllele2(){
        int g = 0;
        int p = beg + 2;
        if (p >= end) 
            return MISSING_GENOTYPE; 
        if (line[p] == '.') 
            return MISSING_GENOTYPE;
        if (line[p] < '0') 
            REPORT("Wrong genotype detected. [2]");
        else 
            g += line[p] - '0';
        return g;
    };
    bool isPhased(){
        if (end - beg != 3) return false;
        int p = beg + 1;
        if (line[p] == '/') return false;
        if (line[p] == '|') return true;
        return false;
    };
    bool isHaploid(){
        return (end - beg == 1);
    };
  private:
    std::string retStr; // this held for the return value;
};

int parseTillChar(const char c, const char* line, const int beg, VCFValue* ret) {
    assert(ret);
    ret->line = line;
    ret->beg = beg;
    ret->end = ret->beg;
    while(line[ret->end] != c && line[ret->end] != '\0') {
        ret->end++;
    }
    return ret->end;
}

int parseTillChar(const char* c, const char* line, const int beg, VCFValue* ret) {
    assert(ret);
    ret->line = line;
    ret->beg = beg;
    ret->end = ret->beg;
    int cLen = strlen(c);
    while(!index(c, line[ret->end])){ // line[beg] is not a separator
        ret->end++;
    }
    return ret->end;
}


class VCFInfoValue{
  public:
    int fingerMark; 
    VCFValue* value;
    VCFInfoValue(){
        this->value = new VCFValue;
    }
    ~VCFInfoValue(){
        delete this->value;
        this->value = NULL;
    }
};

class VCFInfo{
public:
    const char* getTag(const char* tag) {
        if (!tag || tag[0] == '\0') 
            return NULL;

        std::string s = tag;
        if (!hasParsed){
            this->parseActual();
        }
        this->tableIter = this->table.find(s);
        if (this->tableIter != this->table.end()){
            return this->tableIter->second->value->toStr();
        } else {
            return NULL;
        }
    } ;
    ~VCFInfo(){
        for (this->tableIter = this->table.begin(); 
             this->tableIter != this->table.end(); 
             this->tableIter ++){
            if (this->tableIter->second != NULL){
                delete this->tableIter->second;
            }
        }
    };
    void reset() { this-> hasParsed = false;};
    void parse(VCFValue* v) {
        this->value = v;
        this->hasParsed = false;
        this->fingerMark ++ ;
    };
    void parseActual(){
        this->hasParsed = true;
        const char* line = this->value->line;
        int b = this->value->beg;
        int e = this->value->beg;
        static std::string key;

        while ( e < this->value->end){
            key.clear();
            // find tag name;
            while(line[e] != '='){
                key.push_back(this->value->line[e++]);
            }
            b = e + 1; // skip '='
        
            this->tableIter = this->table.find(key);
            if ( this->tableIter == this->table.end()){
                VCFInfoValue* f = new VCFInfoValue;
                this->table[key] = f;
                e = parseTillChar(";\t", line, b, f->value);
                f->fingerMark = this->fingerMark;
            } else {
                e = parseTillChar(";\t", line, b, this->tableIter->second->value);
                this->tableIter->second->fingerMark = this->fingerMark;
            };
            e ++ ;
        }
    };
private:
    bool hasParsed;
    VCFValue* value;
    int fingerMark;
    std::map<std::string, VCFInfoValue*> table;
    std::map<std::string, VCFInfoValue*>::iterator tableIter;
};

typedef int(*PARSE_FUNCTION)(const char* s, int beg, int end);
int parseGenotype(const char* line, int beg, int end){
    if (line[beg] == '.') return  MISSING_GENOTYPE;
    int gt = 0;
    if (!isdigit(line[beg])){
        FATAL("Wrong genotype format");
    };
    gt += (line[beg] - '0');
    if (line[beg+1] == ':')  // sex chromosome
        return gt; 
    if (line[beg+1] == '/' || line[beg+1] == '|'){
        if (line[beg+2] == '.') {
            return MISSING_GENOTYPE;
        }
        if (!isdigit(line[beg+2])){
            FATAL("Wrong genotype format");
        };
        gt += (line[beg+2 ] - '0');
        return gt;
    } else {
        FATAL("Wrong genotype format");
    }
    return -1;
};

// we assume format are always  GT:DP:GQ:GL
class VCFIndividual{
public:

    // FUNC parseFunction[4];
    VCFIndividual(){
        this->include();  // by default, enable everyone
    };
    /**
     * 0-base index for beg and end, e.g.
     *     0 1 2  3
     *     A B C \t
     * beg = 0, end = 3 (line[end] = '\t' or line[end] = '\0')
     */
    int parse(const char* line, const int beg) {
        // need to consider missing genotype
        // need to consider missing field
        
        // skip to next 
        if (!this->isInUse()) {
            return parseTillChar('\t', line, beg, & this->data);
        }
        
        this->data.line = line;
        this->data.beg = beg;
        
        this->fd.clear();
        VCFValue v;
        v.line = line;
        v.beg = beg;
        v.end = beg;
        while (line[v.end] != '\0') { // parse individual values
            v.end = parseTillChar(":\t\0", line, v.beg, &v);
            fd.push_back(v);
            if (line[v.end] == '\t' || line[v.end] == '\0')
                break;
            v.beg = v.end + 1;
        }
        this->data.end = v.end;
        return this->data.end;
    };

    const std::string& getName() const {return this->name;};
    void setName(std::string& s) {this->name = s;};
    void include() {this->inUse = true;};
    void exclude() {this->inUse = false;};
    bool isInUse() {return this->inUse;};
    const VCFValue& operator [] (const unsigned int i) const {
        if (i >= fd.size()){
            FATAL("index out of bound!");
        }
        return (this->fd[i]);                   
    };
    VCFValue& operator [] (const unsigned int i) {
        if (i >= fd.size()){
            FATAL("index out of bound!");
        }
        return (this->fd[i]);
    };
    VCFValue& getData() {return this->data;};
private:
    bool inUse;
    std::string name;

    VCFValue data;            // whole field for the individual
    std::vector<VCFValue> fd; // each field separated by ':'
}; // end VCFIndividual

typedef OrderedMap<int, VCFIndividual*> VCFPeople;

class VCFRecord{
public:
    int parse(const char* line){
        this->line = line;
        // go through VCF sites (first 9 columns)
        int beg = 0;
        int end = 0;
        end = parseTillChar('\t', line, beg, &this->chrom);
        beg = end + 1;
        end = parseTillChar('\t', line, beg, &this->pos);
        beg = end + 1;
        end = parseTillChar('\t', line, beg, &this->id);
        beg = end + 1;
        end = parseTillChar('\t', line, beg, &this->ref);
        beg = end + 1;
        end = parseTillChar('\t', line, beg, &this->alt);
        beg = end + 1;
        end = parseTillChar('\t', line, beg, &this->qual);
        beg = end + 1;
        end = parseTillChar('\t', line, beg, &this->filt);
        beg = end + 1;
        end = parseTillChar('\t', line, beg, &this->info);
        this->vcfInfo.parse(&this->info); // lazy parse inside VCFInfo
        beg = end + 1;
        end = parseTillChar('\t', line, beg, &this->format);
        
        // now comes each individual genotype
        int idx = 0; // peopleIdx
        VCFIndividual* p = this->allIndv[idx]; 
        
        while (line[end] != '\0') {
            beg = end + 1;
            end = p->parse(line, beg);
            if (line[end] == '\0') {
                break;
            }
            idx ++ ;
            if (idx >= this->allIndv.size()){
                FATAL("VCF header and VCF content do not match!");
            }
            p = this->allIndv[idx];
        }
    };
    void createIndividual(const std::string& line){
        std::vector<std::string> sa;
        stringTokenize(line, '\t', &sa);
        if (sa.size() <= 9){
            FATAL("not enough people in the VCF");
        }
        for (int i = 9; i < sa.size(); i++ ) {
            int idx = i - 9;
            VCFIndividual* p = new VCFIndividual;
            this->allIndv[idx] = p;
            p->setName(sa[i]);
        }
    };
    void deleteIndividual(){
        for (int i = 0; i < this->allIndv.size(); i++) {
            if (this->allIndv[i])
                delete this->allIndv[i];
            this->allIndv[i] = NULL;
        }
    };
    void includePeople(const std::string& name){
        for (unsigned int i = 0 ; i != this->allIndv.size() ; i++) {
            VCFIndividual* p = this->allIndv[i];
            if (p->getName() == name) {
                p->include();
            }
        }
    };
    void includePeople(const std::vector<std::string>& v){
        for (unsigned int i = 0; i++ ; i < v.size()){
            this->includePeople(v[i]);
        }
    };
    void includePeopleFromFile(const char* fn){
        LineReader lr(fn);
        std::vector<std::string> fd;
        while(lr.readLineBySep(&fd, "\t ")) {
            for (unsigned int i = 0; i < fd.size(); i++)
                this->includePeople(fd[i]);
        }
    };
    void excludePeople(const std::string& name){
        for (unsigned int i = 0 ; i != this->allIndv.size() ; i++) {
            VCFIndividual* p = this->allIndv[i];
            if (p->getName() == name) {
                p->exclude();
            }
        }
    };
    void excludePeople(const std::vector<std::string>& v){
        for (unsigned int i = 0; i++ ; i != v.size()){
            this->excludePeople(v[i]);
        }
    };
    void excludePeopleFromFile(const char* fn){
        LineReader lr(fn);
        std::vector<std::string> fd;
        while(lr.readLineBySep(&fd, "\t ")) {
            for (unsigned int i = 0; i != fd.size(); i++)
                this->excludePeople(fd[i]);
        }
    };
    const char* getInfoTag(const char* tag) {
        return this->vcfInfo.getTag(tag);
    };
public:
    const char* getChrom() { return this->chrom.toStr(); };
    const int getPos()  { return this->pos.toInt(); };
    const char* getID() { return this->id.toStr(); };
    const char* getRef() { return this->ref.toStr(); };
    const char* getAlt() { return this->alt.toStr(); };
    const char* getQual() { return this->qual.toStr(); };
    const char* getFilt() { return this->filt.toStr(); };
    const char* getInfo() { return this->info.toStr(); };
    const char* getFormat() { return this->format.toStr(); };
    const char* getLine() {return this->line;};
    VCFPeople& getPeople(){
        static bool hasAccess = false;
        if (!hasAccess) {
            for (int i = 0; i < this->allIndv.size(); i++){
                if (allIndv[i]->isInUse()) {
                    this->selectedIndv[this->selectedIndv.size()] = allIndv[i];
                }
            }
            hasAccess = true;
        }
        return this->selectedIndv;
    };
private:
    VCFPeople allIndv;      // all individual
    VCFPeople selectedIndv; // user-selected individual
    VCFValue chrom;
    VCFValue pos;
    VCFValue id;
    VCFValue ref;
    VCFValue alt;
    VCFValue qual;
    VCFValue filt;
    VCFValue info;
    VCFValue format;
    const char* line; // points to data line
    VCFInfo vcfInfo;
}; // VCFRecord

class VCFInputFile{
public:
    VCFInputFile (const char* fn):
        fp(NULL), tabixHandle(NULL){
        this->fileName = fn;
        // open file
        this->fp = new LineReader(fn);
        
        // read header
        while (this->fp->readLine(&this->line)){
            if (line[0] == '#') {
                this->header.push_back(line);
                if (line.substr(0, 6) == "#CHROM") {
                    this->record.createIndividual(line);
                    this->headerLoaded = true;
                    break;
                }
                continue;
            }
            if (line[0] != '#') {
                FATAL("Wrong VCF header");
            }
        }

        this->openIndex();
    };
    ~VCFInputFile(){
        closeIndex();
        this->record.deleteIndividual();
        if (this->fp) delete this->fp;
    };

    void rewriteVCFHeader() {
        std::string s = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        VCFPeople& people = this->record.getPeople();
        for (int i = 0; i <people.size(); i++ ){
            s += '\t';
            s += people[i]->getName();
        }
        this->header[this->header.size()-1] = s;
    };
    VCFHeader* getVCFHeader() {
        this->rewriteVCFHeader();
        return &this->header;
    };
    bool openIndex() {
        return this->openIndex( (fileName + ".tbi").c_str());
    };
    bool openIndex(const char* fn) {
        if (( this->tabixHandle = ti_open(this->fileName.c_str(), 0)) == 0 ) {
            // failed to open tabix index
            this->hasIndex = false;
            return false;
        } else{
			ti_lazy_index_load(this->tabixHandle);
            this->hasIndex = true;
            return true;
        }
    };
    void closeIndex(){
        ti_close(this->tabixHandle);
    };

    void setRangeFile(const char* fn) {
        this->range.addRangeFile(fn);
    }
    void setRangeList(const char* l){
        this->range.addRangeList(l);
    }

    bool readRecord(){
        assert(this->headerLoaded);
        // load contents 
        if (this->range.size() > 0) {
            if (this->hasIndex) {                 // there is index
                static int rangeIdx = 0;
                // unsigned int numRange = this->range->size();
                static ti_iter_t iter; 
                static const char* s = 0;
                int len;
                while (rangeIdx < this->range.size()) {
                    if (!s) { // last time does not read a valid line
                        // get range
                        std::string r;
                        this->range.obtainRange(rangeIdx, &r);
                        // parse range
                        int tid, beg, end;
                        if (ti_parse_region(tabixHandle->idx, r.c_str(), &tid, &beg, &end) != 0){
                            FATAL("Cannot ti_parse_region");
                        }
                        iter =  ti_queryi(tabixHandle, tid, beg, end);
                        s = ti_read(this->tabixHandle, iter, &len);
                        if (s) { // s is valid
                            this->line = s;
                            this->record.parse(this->line.c_str());
                            return true;
                        } else{
                            rangeIdx ++;
                            continue;
                        }
                    } else {  // last time read a valid line
                        s = ti_read(this->tabixHandle, iter, &len);
                        if (!s) {
                            rangeIdx ++;
                            continue;
                        } else {
                            this->line = s;
                            this->record.parse(line.c_str());
                            return true;
                        }
                    }
                } // end while 
                return false;
            } // end hasIndex
            else { // no index
                while(this->fp->readLine(&this->line)){
                    this->record.parse(line.c_str());
                    if (!this->range.isInRange(this->record.getChrom(), this->record.getPos()))
                        continue;
                }
                this->record.parse(this->line.c_str());
            };
        } else { // go line by line
            if (this->fp->readLine(&this->line)){
                this->record.parse(line.c_str());
                return true;
            } else {
                return false; // reach file end
            }
        }
    };
    VCFRecord& getVCFRecord() {return this->record;};
    const std::string& getLine() const {return this->line;};
private:
    VCFHeader header;
    VCFRecord record;
    
    LineReader* fp;
    tabix_t * tabixHandle;
    RangeList range;
    
    std::string fileName;
    std::string line; // actual data line

    bool headerLoaded;
    bool hasIndex;
};

class VCFOutputFile{
public:
    VCFOutputFile(const char* fn) {
        this->fp = new FileWriter(fn);
        if (!this->fp){
            REPORT("Cannot create VCF file!");
            abort();
        }
    };
    ~VCFOutputFile(){
        if (this->fp){
            delete this->fp;
            this->fp = NULL;
        }
    };
    void writeHeader(const VCFHeader* h){
        for (int i = 0; i < h->size(); i++){
            this->fp->writeLine(h->at(i).c_str());
        }
    };
    void writeRecord(VCFRecord* r){
        this->fp->printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
                         r->getChrom(),
                         r->getPos(),
                         r->getID(),
                         r->getRef(),
                         r->getAlt(),
                         r->getQual(),
                         r->getInfo(),
                         r->getFilt(),
                         r->getFormat());
        VCFPeople& p = r->getPeople();
        for (int i = 0; i < p.size() ; i ++ ) {
            VCFIndividual* indv = p[i];
            this->fp->printf("\t%s", indv->getData().toStr());
        }
        this->fp->printf("\n");
    };
private:
    FileWriter* fp;
};

/****************************/
/*    Binary PLINK format   */
/****************************/
/*
 * Documentation 
 * http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml 
 * BED (binary PED):
 * BIM (extended MAP file): chromosome, SNP, cM, base-position, allele 1, allele 2
 *    e.g.
 *    1       snp1    0       1       G       A
 * FAM (first 6 columns of PED file)
 *      Family ID
        Individual ID
        Paternal ID
        Maternal ID
        Sex (1=male; 2=female; other=unknown)
        Phenotype
 *    e.g.
 *    1 1 0 0 1 0
 *    

 */
class PlinkOutputFile{
public:
    PlinkOutputFile(const char* fnPrefix) {
        std::string prefix = fnPrefix;
        this->fpBed = fopen( (prefix + ".bed").c_str(), "wb");
        this->fpBim = fopen( (prefix + ".bim").c_str(), "wt");
        this->fpFam = fopen( (prefix + ".fam").c_str(), "wt");
        if (!this->fpBed || !this->fpBim || !this->fpFam){
            REPORT("Cannot create binary PLINK file!");
            abort();
        }
        // write Bed header
        char c;
        // magic number
        c = 0x6c; // 0b01101100;
        fwrite(&c, sizeof(char), 1, this->fpBed);
        c = 0x1b; // 0b00011011;
        fwrite(&c, sizeof(char), 1, this->fpBed);
        // snp major mode
        c = 0x01; //0b00000001;
        fwrite(&c, sizeof(char), 1, this->fpBed);
    };
    ~PlinkOutputFile() {
        fclose(this->fpBed);
        fclose(this->fpBim);
        fclose(this->fpFam);
    };
    void writeHeader(const VCFHeader* h){
        std::vector<std::string> people; 
        h->getPeopleName(&people);
        this->writeFAM(people);
    };
    // @pos is from 0 to 3
    void setGenotype(unsigned char* c, const int pos, const int geno){
        (*c) |= (geno << (pos<<1));
    }

    void writeRecord(VCFRecord* r){
        // write BIM
        printf("id= %s and its address id = %p\n", r->getID(), r->getID());
        this->writeBIM(r->getChrom(), r->getID(), 0, r->getPos(), r->getRef(), r->getAlt());
        printf("id= %s and its address id = %p\n", r->getID(), r->getID());

        // write BED
        VCFPeople& people = r->getPeople();
        unsigned char c = 0;
        VCFIndividual* indv;
        int offset;
        for (int i = 0; i < people.size() ; i ++) {
            indv = people[i];
            offset = i & (4 - 1);
            if ((*indv)[0].isHaploid()) {
                int a1 = (*indv)[0].getAllele1();
                if (a1 == 0) 
                    setGenotype(&c, offset, HOM_REF);
                else if (a1 == 1)
                    setGenotype(&c, offset, HET);
                else
                    setGenotype(&c, offset, MISSING);
            } else {
                int a1 = (*indv)[0].getAllele1();
                int a2 = (*indv)[0].getAllele2();
                if (a1 == 0) {
                    if (a2 == 0) {
                        //homo ref: 0b00
                    } else if (a2 == 1) {
                        setGenotype(&c, offset, HET); // het: 0b01
                    } else {
                        setGenotype(&c, offset, MISSING); // missing 0b10
                    }
                } else if (a1 == 1) {
                    if (a2 == 0) {
                        setGenotype(&c, offset, HET); // het: 0b01
                    } else if (a2 == 1) {
                        setGenotype(&c, offset, HOM_ALT); // hom alt: 0b11
                    } else {
                        setGenotype(&c, offset, MISSING); // missing
                    }
                }
            }
            if ( offset == 3) { // 3: 4 - 1, so every 4 genotype we will flush 
                fwrite(&c, sizeof(char), 1, this->fpBed);
                c = 0;
            }
        };
        if (offset)
            fwrite(&c, sizeof(char), 1, this->fpBed);
    }
    void writeBIM(const char* chr, const char* id, int mapDist, int pos, const char* ref, const char* alt){
        printf("In writeBIM(), id = %s and its address is id = %p \n", id, id);
        if (strlen(ref) > 1 || strlen(alt) > 1) {
            fprintf(stdout, "skip with ref = %s and alt = %s\n", ref, alt);
            return;
        }
        std::string chrom = chr;
        if (atoi(chr) > 0) {
            fputs(chr, this->fpBim);
            fputc('\t', this->fpBim);
        } else if (chrom == "X")
            fputs("23\t", this->fpBim);
        else if (chrom == "Y")
            fputs("24\t", this->fpBim);
        else if (chrom == "MT")
            fputs("25\t", this->fpBim);
        else {
            fprintf(stdout, "skip chrom %s\n", chr);
            return;
        }
        if (id && id[0] != '.')
            fprintf(this->fpBim, "%s\t", id);
        else
            fprintf(this->fpBim, "%s:%d\t", chrom.c_str(), pos);

        fprintf(this->fpBim, "0\t");
        fprintf(this->fpBim, "%d\t", pos);
        fprintf(this->fpBim, "%c\t", ref[0]);
        fprintf(this->fpBim, "%c\n", alt[0]);
    };
    void writeFAM(std::vector< std::string >& people){
        for (int i = 0; i < people.size(); i++) {
            fprintf(this->fpFam, "%s\t%s\t0\t0\t0\t-9\n", people[i].c_str(), people[i].c_str());
        };
    };
    // NOTE: m should be: marker x people
    void writeBED(Matrix* mat){ 
        int nPeople = mat->cols;
        int nMarker = mat->rows;
        unsigned char c = 0;
        int offset;
        for (int m = 0; m < nMarker; m++){
            for (int i = 0; i < nPeople ; i ++) {
                offset = i & (4 - 1);
                int geno = (int)((*mat)[m][i]);
                switch(geno){
                case HOM_REF:
                    setGenotype(&c, offset, HET); // het: 0b01                
                    break;
                case HET:
                    setGenotype(&c, offset, HET); // het: 0b01
                    break;
                case HOM_ALT:
                    setGenotype(&c, offset, HOM_ALT); // hom alt: 0b11
                    break;
                default:
                    setGenotype(&c, offset, MISSING); // missing
                    break;
                }
            }
            if ( offset == 3) { // 3: 4 - 1, so every 4 genotype we will flush 
                fwrite(&c, sizeof(char), 1, this->fpBed);
                c = 0;
            }
        };
        if (offset)
            fwrite(&c, sizeof(char), 1, this->fpBed);
    };
private:
    // we reverse the two bits as defined in PLINK format, 
    // so we can process 2-bit at a time.
    const static unsigned char HOM_REF = 0x0;     //0b00;
    const static unsigned char HET = 0x2;         //0b10;
    const static unsigned char HOM_ALT = 0x3;     //0b11;
    const static unsigned char MISSING = 0x1;     //0b01;


    FILE* fpBed;
    FILE* fpBim;
    FILE* fpFam;
};


#endif /* _VCFUTIL_H_ */
