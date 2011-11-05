/**
   immediately TODO:
   1. suppport PLINK output
   2. support access INFO tag
   futher TODO:
   1. handle different format GT:GD:DP ...
   2. easy argument processing: clean the argument processing codes.
*/
#include "Argument.h"
#include "IO.h"
#include "./tabix-0.2.2/tabix.h"

#include <cassert>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "PeopleSet.h"
#include "RangeList.h"
#include "Utils.h"

void REPORT(const char* x) { 
    fprintf(stderr, "Report '%s' to zhanxw@umich.edu\n", x ); 
};

void FATAL(const char* x) {
    REPORT(x);
    abort();
};

void REQUIRE_STRING_PARAMETER(const std::string& flag, const char* msg){
    if (flag.size() == 0){
        fprintf(stderr, "%s\n", msg);
        abort();
    }
};

typedef std::vector<std::string> VCFHeader;

#define MISSING_GENOTYPE -1
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
    std::string toStr() const { 
        std::string s;
        for (int i = beg; i < end; i++){
            s.push_back(line[i]);
        }
        return s;
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
            REPORT("Wrong genotype detected. [1]");
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
            REPORT("Wrong genotype detected. [1]");
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
    VCFIndividual():
        isMasked(true)  // by default, enable every one
        {
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
        if (this->isMasked) {
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
            if (line[v.end] == '\t' || line[v.end] == '\0')
                break;
            fd.push_back(v);
            v.beg = v.end + 1;
        }
        this->data.end = v.end;
        return this->data.end;
    };

    const std::string& getName() const {return this->name;};
    void setName(std::string& s) {this->name = s;};
    void addMask() {this->isMasked = true;};
    void delMask() {this->isMasked = false;};
    bool hasMask() {return this->isMasked;};
    const VCFValue& operator [] (const unsigned int i) const {
        return (this->fd[i]);
    };
    VCFValue& getData() {return this->data;};
private:
    bool isMasked;
    std::string name;

    VCFValue data;            // whole field for the individual
    std::vector<VCFValue> fd; // each field separated by ':'
};

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
        end = parseTillChar('\t', line, beg, &this->info);
        beg = end + 1;
        end = parseTillChar('\t', line, beg, &this->filt);
        beg = end + 1;
        end = parseTillChar('\t', line, beg, &this->format);
        
        // now comes each individual genotype
        int idx = 0; // peopleIdx
        VCFIndividual* p = this->allIndv[idx]; 
        
        while (line[end] != '\0') {
            beg = end + 1;
            end = p->parse(line, beg);
            if (line[end] == '\0') 
                break;
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
            delete this->allIndv[i];
            this->allIndv[i] = NULL;
        }
    };
    void setPeople(PeopleSet* pInclude, PeopleSet* pExclude) {
        if (pInclude && pInclude->size() > 0) {
            // by default set the mask, unless the people is in pInclude and not in pExclude
            for (unsigned int i = 0; i != this->allIndv.size(); i++){
                VCFIndividual* p = this->allIndv[i];
                p->addMask();
                if (pInclude->contain(p->getName())) {
                    if (!pExclude)
                        p->delMask();
                    else if (!pExclude->contain(p->getName())) {
                        p->delMask();
                    }
                }
            }
        } else {
            // by default clear the mask, unless the people is in pExclude
            for (unsigned int i = 0; i != this->allIndv.size(); i++){
                VCFIndividual* p = this->allIndv[i];
                p->delMask();
                if (!pExclude && pExclude->contain(p->getName()))
                    p->addMask();
            }
        }
    };
public:
    const std::string getChrom() const { return this->chrom.toStr(); };
    const int getPos() const { return this->pos.toInt(); };
    const std::string getID() const { return this->id.toStr(); };
    const std::string getRef() const { return this->ref.toStr(); };
    const std::string getAlt() const { return this->alt.toStr(); };
    const std::string getQual() const { return this->qual.toStr(); };
    const std::string getFilt() const { return this->filt.toStr(); };
    const std::string getInfo() const { return this->info.toStr(); };
    const std::string getFormat() const { return this->format.toStr(); };
    const char* getLine() const {return this->line;};
    VCFPeople& getPeople(){
        static bool hasAccess = false;
        if (!hasAccess) {
            for (int i = 0; i < this->allIndv.size(); i++){
                if (!allIndv[i]->hasMask()) {
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
}; // VCFRecord

class VCFInputFile{
public:
    VCFInputFile (const char* fn):
        fp(NULL), tabixHandle(NULL), range(NULL){
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
            this->hasIndex = false;
            return false;
        } else{
            this->hasIndex = true;
            return true;
        }
    };
    void closeIndex(){
        ti_close(this->tabixHandle);
    };

    void setRange(RangeList* rl) { this->range = rl; };
    void setPeople(PeopleSet* pInclude, PeopleSet* pExclude) {
        this->record.setPeople(pInclude, pExclude);
    };
    bool readRecord(){
        assert(this->headerLoaded);
        // load contents 
        if (this->range && this->range->size() > 0) {
            if (this->hasIndex) {                 // there is index
                static int rangeIdx = 0;
                // unsigned int numRange = this->range->size();
                static ti_iter_t iter; 
                static const char* s = 0;
                int len;
                while (rangeIdx < this->range->size()) {
                    if (!s) { // last time does not read a valid line
                        // get range
                        std::string r;
                        this->range->obtainRange(rangeIdx, &r);
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
                            continue;
                        }
                    } else {  // last time read a valid line
                        s = ti_read(this->tabixHandle, iter, &len);
                        if (!s) {
                            rangeIdx ++;
                            continue;
                        }
                    }
                } // end while 
                return true;
            } // end hasIndex
            else { // no index
                while(this->fp->readLine(&this->line)){
                    this->record.parse(line.c_str());
                    if (!this->range->isInRange(this->record.getChrom(), this->record.getPos()))
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
    RangeList* range;
    
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
    void writeHeader(const VCFHeader* h){
        for (int i = 0; i < h->size(); i++){
            this->fp->writeLine(h->at(i).c_str());
        }
    };
    void writeRecord(VCFRecord* r){
        this->fp->printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
                         r->getChrom().c_str(),
                         r->getPos(),
                         r->getID().c_str(),
                         r->getRef().c_str(),
                         r->getAlt().c_str(),
                         r->getQual().c_str(),
                         r->getInfo().c_str(),
                         r->getFilt().c_str(),
                         r->getFormat().c_str());
        VCFPeople& p = r->getPeople();
        for (int i = 0; i < p.size() ; i ++ ) {
            VCFIndividual* indv = p[i];
            this->fp->printf("\t%s", indv->getData().toStr().c_str());
        }
        this->fp->printf("\n");
    };
    ~VCFOutputFile(){
        if (this->fp){
            delete this->fp;
        }
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
    };
    ~PlinkOutputFile() {
        fclose(this->fpBed);
        fclose(this->fpBim);
        fclose(this->fpFam);
    };
    void writeHeader(const VCFHeader* h){
        std::string lastLine = h->at(h->size() - 1);
        std::vector<std::string> people; 
        stringTokenize(lastLine, "\t", &people);
        if (people.size() < 9) {
            fprintf(stderr, "Wrong VCF Header line: %s\n", lastLine.c_str());
        };
        for (int i = 9; i < people.size(); i++) {
            fprintf(this->fpFam, "%s\t\"");
        };
    };
    void writeRecord(VCFRecord* r){
private:
    FILE* fpBed;
    FILE* fpBim;
    FILE* fpFam;
};

int main(int argc, char** argv){
    time_t currentTime = time(0);
    fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

    ////////////////////////////////////////////////
    BEGIN_PARAMETER_LIST(pl)
        ADD_PARAMETER_GROUP(pl, "Input/Output")
        ADD_STRING_PARAMETER(pl, input, "--input", "input VCF File")
        ADD_STRING_PARAMETER(pl, output, "--output", "output prefix")
        ADD_PARAMETER_GROUP(pl, "People Filter")
        ADD_STRING_PARAMETER(pl, peopleIncludeID, "--peopleIncludeID", "give IDs of people that will be included in study")
        ADD_STRING_PARAMETER(pl, peopleIncludeFile, "--peopleIncludeFile", "from given file, set IDs of people that will be included in study")
        END_PARAMETER_LIST(pl)
        ;    

    pl.Read(argc, argv);
    pl.Status();

    REQUIRE_STRING_PARAMETER(FLAG_input, "Please provide input file using: --input");
    REQUIRE_STRING_PARAMETER(FLAG_output, "Please provide output prefix using: --output");

    const char* fn = FLAG_input.c_str(); //"test.vcf";
    VCFInputFile vin(fn);

    // set range filters here
    RangeList rl;

    // set people filters here
    PeopleSet peopleInclude;
    PeopleSet peopleExclude;
    peopleInclude.readID(FLAG_peopleIncludeID.c_str());
    vin.setPeople(&peopleInclude, &peopleExclude);

    // let's write it out.
    const char* fout = "test.out.vcf";
    VCFOutputFile vout(fout);
    vout.writeHeader(vin.getVCFHeader());
    while (vin.readRecord()){
        VCFRecord& r = vin.getVCFRecord(); 
        VCFPeople& people = r.getPeople();
        const VCFIndividual* indv;
        vout.writeRecord(& r);
//        printf("%s:%d\t", r.getChrom().c_str(), r.getPos());
//        for (int i = 0; i < people.size(); i++) {
//            indv = people[i];
//            printf("%d ", (*indv)[0].toInt());  // [0] meaning the first field of each individual
//        }
//        printf("\n");
    };

    currentTime = time(0);
    fprintf(stderr, "Analysis ended at: %s", ctime(&currentTime));

    return 0;
};
