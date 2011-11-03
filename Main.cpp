/**
   TODO:
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

#include "InfoGrepper.h"
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
typedef enum {VCF_CHROM= 0,
              VCF_POS  = 1,
              VCF_ID   = 2,
              VCF_REF  = 3,
              VCF_ALT  = 4,
              VCF_QUAL = 5,
              VCF_FILTER = 6,
              VCF_INFO   = 7,
              VCF_FORMAT = 8,
              VCF_UNDEF = -1} 
    VCF_SITE_t;

typedef std::vector<std::string> VCFHeader;

class VCFValue{
public:
    int beg;
    int end;
    const char* line;
public:
    int toInt();
    void toInt(int* i);
    double toDouble();
    void toDouble(double* d);
    std::string toStr();
    void toStr(std::string* s) {};
public:
    int getGenotype(){
    };
    int getAllele1(){
    };
    int getAllele2(){
    };
    bool isPhased(){
    };
};

#define MISSING_GENOTYPE -1
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

int parseDepth(const char* line, int beg, int end){
    char s[32];
    int p = beg;
    while (isdigit(line[p]) && p <= end){
        p++;
    }
    if (p == end)
        FATAL("Wrong depth format");
    strncpy(s, &line[beg], p - beg + 1);
    return atoi(s);
};

int parseGenotypeQuality(const char* line, int beg, int end){
    char s[32];
    int p = beg;
    while (isdigit(line[p]) && p <= end){
        p++;
    }
    if (p == end)
        FATAL("Wrong genotype quality format");
    strncpy(s, &line[beg], p - beg + 1);
    return atoi(s);
};

// we assume format are always  GT:DP:GQ:GL
class VCFIndividual{
public:

    // FUNC parseFunction[4];
    VCFIndividual():
        isMasked(true)  // by default, enable every one
        {
            this->parseFunction[0] = parseGenotype; //parseGenotype;
            this->parseFunction[1] = parseDepth; // TODO: change to other function
            // this->parseFunction[2] = parseGenotype;
            // this->parseFunction[3] = parseGenotype;
        };
    /**
     * 0-base index for beg and end, e.g.
     *     0 1 2  3
     *     A B C \t
     * beg = 0, end = 3 (line[end] = '\t' or line[end] = '\0')
     */
    bool parse(const char* line, int beg, int end) {
        // need to consider missing genotype
        // need to consider missing field
        if (this->isMasked) return true;

        this->beg = beg;
        this->end = end;
        int pt = beg;
        int parseFunctionIdx = 0;
        while(line[pt] != ':' && pt <= end){
            pt ++;
            if (line[pt] == ':' || pt == end){
                switch (parseFunctionIdx){
                case 0:
                    this->gt = (*(parseFunction[parseFunctionIdx]))(line, beg, end);
                    break;
                case 1:
                    this->dp = (*(parseFunction[parseFunctionIdx]))(line, beg, end);
                    break;
                case 2:
                    this->gq = (*(parseFunction[parseFunctionIdx]))(line, beg, end);
                    break;
                default:
                    break;
                }
                parseFunctionIdx ++;
            }
        }
        return true;
    };
    const std::string& getName() const {return this->name;};
    void setName(std::string& s) {this->name = s;};
    void addMask() {this->isMasked = true;};
    void delMask() {this->isMasked = false;};
    void resetContent() { gt = -1; dp = -1; gq = -1; };
    int getGenotype() const {return this->gt;};
private:
    bool isMasked;
    std::string name;
    int beg;
    int end;
    PARSE_FUNCTION parseFunction[4];
    int gt;
    int dp;
    int gq;
    
    // int getTagStr(const char* tag, char* value, int len);
    // int getTagInt(const char* tag, int* value);
    // int getTagDouble(const char* tag, double* value);
    // int getGenotype(const char* tag, int* geno);
};

typedef OrderedMap<int, VCFIndividual*> VCFPeople;

class VCFRecord{
public:
    int parse(const char* line){
        this->line = line;
        // go through each character
        int beg = 0;
        int end = 0;
        VCF_SITE_t state = VCF_CHROM;
        while (line[end] != '\0'){
            end ++ ;
            if (line[end] == '\t'){
                // update states accordingly
                std::string s;
                switch (state){
                case VCF_CHROM:
                    this->chrom.assign(& line[beg], &line[end]);
                    state = VCF_POS;
                    break;
                case VCF_POS:
                    s.assign(& line[beg], &line[end]);
                    if (!str2int(s.c_str(), &this->pos)) {
                        FATAL("Position error!");
                    }
                    state = VCF_ID;
                    break;
                case VCF_ID:
                    this->id.assign(& line[beg], &line[end]);
                    state = VCF_REF;
                    break;
                case VCF_REF:
                    this->ref.assign(& line[beg], &line[end]);
                    state = VCF_ALT;
                    break;
                case VCF_ALT:
                    this->alt.assign(& line[beg], &line[end]);
                    state = VCF_QUAL;
                    break;
                case VCF_QUAL:
                    this->qual.assign(& line[beg], &line[end]);
                    state = VCF_FILTER;
                    break;
                case VCF_FILTER:
                    this->filt.assign(& line[beg], &line[end]);
                    state = VCF_INFO;
                    break;
                case VCF_INFO:
                    this->info.assign(& line[beg], &line[end]);
                    state = VCF_FORMAT;
                    break;
                case VCF_FORMAT:
                    this->format.assign(& line[beg], &line[end]);
                    state = VCF_UNDEF;
                    break;
                default:
                    REPORT("Should not reach here!");
                    break;
                };
                
                beg = end + 1;
                end = beg;
                if (state == VCF_UNDEF) 
                    break;
            }
        }
        // now comes each individual genotype
        int idx = 0; // peopleIdx
        VCFIndividual* p = this->indv[idx]; 
        p->resetContent();

        while (line[end] != '\0') {
            end ++ ;
            if (line[end] == '\t' || line[end] == '\0') {
                if (!p->parse(line, beg, end)){
                    REPORT("Parse indv failed!");
                    abort();
                }
                if (line[end] == '\0') 
                    break;
                beg = end + 1;
                end = beg;
                idx ++ ;
                if (idx >= this->indv.size()){
                    FATAL("VCF header and VCF content do not match!");
                }
                p = this->indv[idx];
                p->resetContent();
            }
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
            this->indv[idx] = p;
            p->setName(sa[i]);
        }
    };
    void deleteIndividual(){
        for (int i = 0; i < this->indv.size(); i++) {
            delete this->indv[i];
            this->indv[i] = NULL;
        }
    };
    void setPeople(PeopleSet* pInclude, PeopleSet* pExclude) {
        if (pInclude && pInclude->size() > 0) {
            // by default set the mask, unless the people is in pInclude and not in pExclude
            for (unsigned int i = 0; i != this->indv.size(); i++){
                VCFIndividual* p = this->indv[i];
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
            for (unsigned int i = 0; i != this->indv.size(); i++){
                VCFIndividual* p = this->indv[i];
                p->delMask();
                if (!pExclude && pExclude->contain(p->getName()))
                    p->addMask();
            }
        }
    };
public:
    const std::string& getChrom() const { return this->chrom; };
    const int getPos() const { return this->pos; };
    const std::string& getID() const { return this->id; };
    const std::string& getRef() const { return this->ref; };
    const std::string& getAlt() const { return this->alt; };
    const std::string& getQual() const { return this->qual; };
    const std::string& getFilt() const { return this->filt; };
    const std::string& getInfo() const { return this->info; };
    const std::string& getFormat() const { return this->format; };
    const char* getLine() const {return this->line;};
    const VCFPeople& getIndv() const { return this->indv;};
private:
    VCFPeople indv; // each individual
    std::string chrom;
    int pos;
    std::string id;
    std::string ref;
    std::string alt;
    std::string qual;
    std::string filt;
    std::string info;
    std::string format;
    const char* line;
};

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

    VCFHeader* getHeader() {return &this->header;};
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
        // if (!this->headerLoaded) {
        //     // load header
        //     this->fp->readByLine(&this->line);
            
        //     if (this->range) {
        //         // check if there is index
        //         if (!this->hasIndex){
        //             fprintf(stderr, "creating index is recommended\n");
        //         }
        //     }
        // } else{
        // load contents 
        if (this->range) {
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
    const VCFRecord& getVCFRecord() const {return this->record;};
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
        this->fp = fopen(fn, "wt");
        if (!this->fp){
            REPORT("Cannot create VCF file!");
            abort();
        }
    };
    void writeHeader(const VCFHeader* h){
        for (int i = 0; i < h->size(); i++){
            fprintf(this->fp, "%s\n", h->at(i).c_str());
        }
    };
    void writeRecord(const VCFRecord* r){
        fprintf(this->fp, "%s\n", r->getLine());
    };
    ~VCFOutputFile(){
        if (this->fp){
            fclose(this->fp);
        }
    };
private:
    FILE* fp;
};

int main(int argc, char** argv){
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
    
    if (FLAG_input.size() == 0) {
        fprintf(stderr, "use debug parameters\n");
        const char* fn = "test.vcf";
        VCFInputFile vin(fn);
        const char* fout = "test.out.vcf";
        VCFOutputFile vout(fout);
        vout.writeHeader(vin.getHeader());
        
        RangeList rl;
        PeopleSet peopleInclude;
        PeopleSet peopleExclude;
        if (FLAG_peopleIncludeID.size() > 0) {
            peopleInclude.readID(FLAG_peopleIncludeID.c_str());
        }
        
        vin.setPeople(&peopleInclude, &peopleExclude);

        while (vin.readRecord()){
            const VCFRecord& r = vin.getVCFRecord(); 
            const VCFPeople& indv = r.getIndv();
            int idx;
            const VCFIndividual* people;
            vout.writeRecord(& r);
            printf("%s:%d\t", r.getChrom().c_str(), r.getPos());
            for (int i = 0; i < indv.size(); i++) {
                people = indv[i];
                printf("%d ", people->getGenotype());
            }
            printf("\n");
        };
    } else{
        if (FLAG_output.size() == 0){
            fprintf(stderr, "Please provide output prefix\n");
        }
    }
    return 0;
};

// struct VCFFilter {
//     // PeopleIndex* peopleIndex;
//     // RangeList* rangeList;
//     int siteDepthMin;
//     int siteDepthMax;
//     int indvDepthMin;
//     int indvDepthMax;
//     int indvQualMin;
//     InfoGrepper* infoGrepper;
//     VCFFilterArgument():
//         peopleIndex(0),
//         rangeList(0),
//         siteDepthMin(0),
//         siteDepthMax(0),
//         indvDepthMin(0),
//         indvDepthMax(0),
//         indvQualMin(0),
//         infoGrepper(0)
//         {};
// };

// class VCFFile{
// public:
//     int  open(const char* fileName);
//     void output012File(const std::string argOutputPrefix, const VCFFilterArgument& filterArgument);
//     void close();
// private:
//     void vcfLineTo012(const char* s, const VCFFilterArgument& filterArgument, std::vector<short int>* geno, FILE* genoFile, FILE* posFile, FILE* frqFile);
//     int  loopEachLine(const VCFFilterArgument& filterArgument, FILE* genoFile, FILE* posFile, FILE* frqFile);
//     void filterByPeople(PeopleIndex& pi) {
//         pi.readVCFheader(this->vcfHeader[this->vcfHeader.size()-1].c_str());
//         pi.setFilteredID(& this->inclusionColumn);
//         // debug code
//         // printf("filtered result: \n");
//         // for (unsigned int i = 0; i < inclusionColumn.size(); i++) {
//         //     printf("index[%d], col = %d\n", i, inclusionColumn[i]);
//         // }
//     };
//     bool indexFileExist() { return true;};
//     // open sub-indexed files
//     // .people  format: col peopleID
//     // .af      format: TODO(zhanxw)
//     int createIndexFiles(const char* fileName);
//     int openIndexFiles(const char* fileName);

// private:
//     IFILE vcfHandle;
//     tabix_t * tabixHandle;
//     // std::vector<bool>    inclusionColumn; // every column in the VCF corresponds to a bool value, only col > 9 has meanings.
//     VCFHeader* header;
//     // IFILE vcfGTHandle;  // store GT
//     // IFILE vcfGQHandle;  // store GP
//     // IFILE vcfDPHandle;  // store DP
// };

// void VCFFile::output012File(const std::string argOutputPrefix, const VCFFilterArgument& filterArgument){
//     // dump indv file
//     filterByPeople(*filterArgument.peopleIndex);
//     std::string fn = argOutputPrefix + ".012.indv";
//     FILE* indvFile = fopen(fn.c_str(), "w");
//     filterArgument.peopleIndex->dumpIndvFile(this->vcfHeader[this->vcfHeader.size()-1].c_str(), indvFile);
//     fclose(indvFile);

//     // open freq, pos file
//     fn = argOutputPrefix + ".012.frq";
//     FILE* freqFile = fopen(fn.c_str(), "w");
//     fn = argOutputPrefix + ".012.pos";
//     FILE* posFile = fopen(fn.c_str(), "w");
//     fn = argOutputPrefix + ".012";
//     FILE* genoFile = fopen(fn.c_str(), "w");

//     loopEachLine(filterArgument, genoFile, posFile, freqFile);

//     fclose(freqFile);
//     fclose(posFile);
//     fclose(genoFile);
// }

// int  VCFFile::loopEachLine(const VCFFilterArgument& filterArgument, FILE* genoFile, FILE* posFile, FILE* frqFile){
//     // parse vcf file content
//     //#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ID1 ...
//     //
//     int tid, beg, end;
//     //TODO: make this flexible
//     unsigned int numRange = filterArgument.rangeList->size();
//     std::string range;

//     const char *s;
//     int len;
//     std::vector<short int> g; // a very large 1-d array store genotype

//     if (numRange == 0) {
//         range = ".";
//         ti_iter_t iter;
//         iter = ti_query(tabixHandle, 0, 0, 0);
//         while ((s = ti_read(tabixHandle, iter, &len)) != 0) {
//             // s is the string we will work on
//             // fputs(s, stdout); fputc('\n', stdout);
//             if ( s[0] == '#') continue;
//             vcfLineTo012(s, filterArgument, &g, genoFile, posFile, frqFile);
//         }
//         ti_iter_destroy(iter);

//     } else {
//         for (unsigned int i = 0; i < numRange; i++) {
//             filterArgument.rangeList->obtainRange(i, &range);
//             if (ti_parse_region(tabixHandle->idx, range.c_str(), &tid, &beg, &end) == 0) {
//                 ti_iter_t iter;
//                 // const char *s;
//                 // int len;
//                 iter = ti_queryi(tabixHandle, tid, beg, end);
//                 while ((s = ti_read(tabixHandle, iter, &len)) != 0) {
//                     // s is the string we will work on
//                     // fputs(s, stdout); fputc('\n', stdout);
//                     vcfLineTo012(s, filterArgument, &g, genoFile, posFile, frqFile);
//                 }
//                 ti_iter_destroy(iter);
//             }
//             else fprintf(stderr, "[main] invalid region: unknown target name or minus interval.\n");
//         }
//     }
//     // transpose genotype matrix
//     long int n = g.size();
//     if (n == 0)
//         return 0;
//     long int nc = 0;
//     //printf("people num (unfilt) = %d\n", this->inclusionColumn.size());
//     for (unsigned int i = 0; i< this->inclusionColumn.size(); i++) {
//         if (this->inclusionColumn[i])
//             ++nc;
//     }
//     long int nr = n / nc;
//     assert( nr * nc == n);

//     if (nc != 1) {
//         transposeMatrix(&g, nr, nc);
//     }
    
//     // output to 012 file
//     std::swap(nr, nc);
//     for (long int i = 0 ; i < nr; i++) {
//         for (long int j = 0; j < nc; j++ ){
//             if (j) fputc('\t', genoFile);
//             if (g[i*nc+j]>=0)
//                 fputc('0' + g[i*nc+j], genoFile);
//             else 
//                 fputs("-1", genoFile);
//         }
//         fputc('\n', genoFile);
//     }
//     return 0;
// }

// void VCFFile::vcfLineTo012(const char* s, const VCFFilterArgument& filterArgument, std::vector<short int>* geno, FILE* genoFile, FILE* posFile, FILE* frqFile){
//     unsigned int state = 0;
//     unsigned int i = 0;
//     unsigned int begin = 0;
//     unsigned int end = 0;

//     int siteDepth = 0;

//     int chrBegin = 0;
//     int posEnd = 0;
//     int idBegin, idEnd;
//     int refBegin, refEnd;
//     int altBegin, altEnd;
//     int freq[2] ={0,0};
//     int all;
//     double f0, f1;
//     int indvQual;
//     int indvDepth;
//     while (s[i++]) {
//         if (s[i] == '\t' || s[i] == '\0') {
//             end = i - 1;
//             // deal with the value
//             switch(state){
//             case VCF_CHROM:
//                 chrBegin = begin;
//                 break;
//             case VCF_POS:
//                 posEnd = end;
//                 break;
//             case VCF_ID:
//                 idBegin = begin;
//                 idEnd = end;
//                 break;
//             case VCF_REF:
//                 refBegin = begin;
//                 refEnd = end;
//                 break;
//             case VCF_ALT:
//                 altBegin = begin;
//                 altEnd = end;
//                 break;
//             case VCF_QUAL:
//                 break;
//             case VCF_FILTER:
//                 break;
//             case VCF_INFO:
//                 siteDepth = obtainDepth(s, begin, end);
//                 if (filterArgument.siteDepthMin != 0 && siteDepth < filterArgument.siteDepthMin){
//                     //printf("depth %d < %d failed\n", depth, argDepth);
//                     goto NEXT_LINE;
//                 }
//                 if (filterArgument.siteDepthMax != 0 && siteDepth > filterArgument.siteDepthMax){
//                     goto NEXT_LINE;
//                 }
//                 if (!filterArgument.infoGrepper->match(s, begin, end)){
//                     //printf("depth %d pass\n", depth);
//                     goto NEXT_LINE;
//                 }
//                 if (!filterArgument.infoGrepper->match(s, begin, end)) {
//                     goto NEXT_LINE;
//                 }
//                 break;
//             case VCF_FORMAT:
//                 break;
//             default:
//                 // check people filter
//                 //printf("column %d is %s", state-9, (this->inclusionColumn[state-9] ? "true": "false"));
//                 if (!this->inclusionColumn[state-9]){
//                     goto NEXT_FIELD;
//                 }

//                 // deal with filtered individuals
//                 //NOTE: only work with GT:DP:GQ:GL format
                
//                 // check GT
//                 //also GT should be X/X format
// #define MISSING_GT_VALUE -1
//                 int gt = 0; // 0, 1, or 2
//                 int g;
//                 if (s[begin] == '.') {
//                     gt = -1;
//                     //fputs(MISSING_GT_VALUE, stdout);
//                     //geno->push_back(-1);
//                 } else {
//                     int temp = s[begin] - '0';
//                     if (temp < 0 || temp > 2) {
//                         REPORT("Genotype out of bound") ;
//                     } else {
//                         g = s[begin] - '0';
//                         gt += g;
//                         freq[g] ++;
//                         // check if it's sex chrom
//                         if (s[begin+1] != ':') {
//                             g = s[begin+2] - '0';
//                             gt += g;
//                             freq[g] ++;
//                         }
//                     }
//                     //fputc('0' + gt, stdout);
//                     // geno->push_back(gt);
//                     //printf("%d", gt);
//                 }
                
//                 // check DP
//                 int innerPos;
//                 if (s[begin+1] != ':') {
//                     innerPos = begin + 3; // innerPos used for position within this field
//                 } else {
//                     innerPos = begin + 1; // innerPos used for position within this field
//                 }
//                 if (s[innerPos] != ':') {
//                     REPORT( "We cannot handle this format");
//                     goto NEXT_LINE;
//                 } else{
//                     int indvDepth = 0;
//                     int depthPos = begin+4;
//                     while (s[depthPos] != ':') {
//                         indvDepth *= 10;
//                         indvDepth += (s[depthPos] - '0');
//                         depthPos++;
//                     }
//                     innerPos = depthPos;
//                     // fprintf(stdout, "indv depth = %d\n", indvDepth);
//                     // filter on indvDepth
//                     if (filterArgument.indvDepthMin != 0 && indvDepth < filterArgument.indvDepthMin){
//                         gt = MISSING_GT_VALUE;
//                     }
//                     if (filterArgument.indvDepthMax != 0 && indvDepth > filterArgument.indvDepthMax){
//                         gt = MISSING_GT_VALUE;
//                     }
//                 }

//                 // check GQ
//                 if (s[innerPos] != ':') {
//                     REPORT( "We cannot handle this format");
//                     goto NEXT_LINE;
//                 } else {                    
//                     int indvQual = 0;
//                     int qualPos = innerPos + 1;
//                     while (s[qualPos] != ':') {
//                         indvQual *= 10;
//                         indvQual += (s[qualPos] - '0');
//                         qualPos++;
//                     }
//                     innerPos = qualPos;
//                     // fprintf(stdout, "indv depth = %d\n", indvDepth);
//                     // filter on indvDepth
//                     if (filterArgument.indvQualMin != 0 && indvQual < filterArgument.indvQualMin){
//                         gt = MISSING_GT_VALUE;
//                     }
//                 }
//                 geno->push_back(gt);
//                 break;
//             } // end switch

//             // finish dealing with the line
//           NEXT_FIELD:
//             ++ state;
//             begin = i+1;
//         }
//     }
//     //fputs("\n", stdout);

// // this macro will print s[b..e] to file f. (boundary inclusive)
// #define FIELD_PRINT(s, b, e, f)                 \
//     do {                                        \
//         for (unsigned int i = b; i <= e; i++)   \
//             fputc(s[i], f);                     \
//     }while(0);                                  \

//     // output pos file
//     FIELD_PRINT(s, chrBegin, posEnd, posFile);
//     fputc('\n', posFile);

//     // output freq file
//     all = freq[0] + freq[1];
//     f0 = 0.0;
//     f1 = 0.0;
//     if (all != 0) {
//         f0 = (double)(freq[0]) / double(all);
//         f1 = 1 - f0;
//     }
//     FIELD_PRINT(s, chrBegin, posEnd, frqFile);
//     fputc('\t', frqFile);
//     FIELD_PRINT(s, idBegin, idEnd, frqFile);
//     fputc('\t', frqFile);
//     FIELD_PRINT(s, refBegin, refEnd, frqFile);
//     fputc(':', frqFile);
//     fprintf(frqFile, "%lf", f0);
//     fputc('\t', frqFile);
//     FIELD_PRINT(s, altBegin, altEnd, frqFile);
//     fputc(':', frqFile);
//     fprintf(frqFile, "%lf", f1);
//     fputc('\n', frqFile);
//   NEXT_LINE:
//     return;
// };

// // return 0 for success
// int VCFFile::open(const char* fileName) {
//     this->vcfHandle = ifopen(fileName, "r");
//     if (!this->vcfHandle) return -1;
//     // read header
//     String line;
//     uint32_t numField = 0;
//     uint32_t numLine = 0;
//     while (line.ReadLine(this->vcfHandle) >= 0) {
//         ++ numLine ;
//         // vcf header part
//         if (line[0] == '#') {
//             this->vcfHeader.push_back(line.c_str());
//         } else {
//             break;
//         }
//     }

//     // check whether indexing files exist
//     if (!indexFileExist()) {
//         // build index
//     }

//     // open index
//     if (( this->tabixHandle = ti_open(fileName, 0)) == 0 ) {
//         fprintf(stderr, "Cannot open use tabix to open: %s\n", fileName);        
//         exit(1);
//     }
//     if (ti_lazy_index_load(tabixHandle)){
//         fprintf(stderr, "Cannot open tabix index file\n");
//         fprintf(stderr, "Use this command to create the tabix index file:\n");
//         fprintf(stderr, "Use this command to create the tabix index file:\n");
//         fprintf(stderr, "(grep ^\"#\" %s; grep -v ^\"#\" %s | sort -k1,1 -k2,2n) | bgzip > sorted.%s.gz;\n", fileName, fileName, fileName);
//         fprintf(stderr, "tabix -p vcf sorted.%s.gz;\n", fileName);
//         exit(1);
//     };

//     //openIndexFiles(fileName);
//     return 0;
// }

// void VCFFile::close() {
//     if (this->vcfHandle)
//         ifclose(this->vcfHandle);

//     // close index
//     ti_close(this->tabixHandle);
// }

// int main(int argc, char *argv[])
// {
//     // parse argumetn
//     // --peopleIncludeID
//     // --peopleIncludeFile
//     // --peopleExcludeID
//     // --peopleExcludeFile
//     // --outputPrefix
//     // --outputFormat [plink, 012]
//     // --geneIncludeName
//     // --geneIncludeFile
//     // --geneExcludeName
//     // --geneExcludeFile
//     // to add
//     // filter by AF,
//     // filter by QUAL
//     // filter by Flanking sequences?

//     String argInputVCFFileName;
//     String argOutputPrefix;

//     String argPeopleIncludeID;
//     String argPeopleIncludeFile;
//     String argPeopleExcludeID;
//     String argPeopleExcludeFile;

//     String argGeneListFile;
//     String argGeneIncludeName;
//     String argRangeList;
//     String argRangeFile;
//     String argSiteDepthMin;
//     String argSiteDepthMax;

//     String argIndvDepthMin;
//     String argIndvDepthMax;
//     String argIndvQualMin;

//     String argInfoGrep;

//     ParameterList pl;
//     BEGIN_LONG_PARAMETERS(longParameters)
//         LONG_PARAMETER_GROUP("Input/Output")
//         LONG_STRINGPARAMETER("input",&argInputVCFFileName)
//         LONG_STRINGPARAMETER("outputPrefix",&argOutputPrefix)

//         LONG_PARAMETER_GROUP("People Filter")
//         LONG_STRINGPARAMETER("peopleIncludeID",&argPeopleIncludeID)
//         LONG_STRINGPARAMETER("peopleIncludeFile",&argPeopleIncludeFile)
//         LONG_STRINGPARAMETER("peopleExcludeID",&argPeopleExcludeID)
//         LONG_STRINGPARAMETER("peopleExcludeFile",&argPeopleExcludeFile)

//         LONG_PARAMETER_GROUP("Site Filter")
//         LONG_STRINGPARAMETER("geneListFile", &argGeneListFile)
//         LONG_STRINGPARAMETER("geneIncludeName", &argGeneIncludeName)
//         LONG_STRINGPARAMETER("rangeList", &argRangeList)
//         LONG_STRINGPARAMETER("rangeFile", &argRangeFile)
//         LONG_STRINGPARAMETER("siteDepthMin", &argSiteDepthMin)
//         LONG_STRINGPARAMETER("siteDepthMax", &argSiteDepthMax)

//         LONG_PARAMETER_GROUP("Individual Filter")
//         LONG_STRINGPARAMETER("indvDepthMin", &argIndvDepthMin)
//         LONG_STRINGPARAMETER("indvDepthMax", &argIndvDepthMax)
//         LONG_STRINGPARAMETER("indvQualMin", &argIndvQualMin)

//         LONG_PARAMETER_GROUP("INFO field Grepper")
//         LONG_STRINGPARAMETER("infoGrep", &argInfoGrep)
//         END_LONG_PARAMETERS();

//     pl.Add(new LongParameters("Available Options", longParameters));
//     pl.Read(argc, argv);
//     pl.Status();

//     if (argInputVCFFileName.Length() == 0) {
//         fprintf(stderr, "Please specify --input\n");
//         exit(1);
//     }
//     if (argOutputPrefix.Length() == 0) {
//         fprintf(stderr, "Please specify --outputPrefix\n");
//         exit(1);
//     }

//     PeopleIndex pi;
//     pi.includeID(argPeopleIncludeID);
//     pi.includeFile(argPeopleIncludeFile);
//     pi.excludeID(argPeopleExcludeID);
//     pi.excludeFile(argPeopleExcludeFile);

//     RangeList rl;
//     if (argGeneIncludeName.Length() !=0 && argGeneListFile.Length() == 0){
//         fprintf(stderr, "Please specify geneListFile, or we don't know where a gene starts and ends. \n");
//     }
//     rl.filterGeneName(argGeneIncludeName, argGeneListFile);
//     if (argRangeList.Length() != 0)
//         rl.addRangeList(argRangeList);
//     if (argRangeFile.Length() != 0)
//         rl.addRangeFile(argRangeFile);

//     InfoGrepper ig;
//     if (argInfoGrep.Length() > 0)
//         ig.readPatter(argInfoGrep);

//     VCFFilterArgument filterArgument;
//     filterArgument.peopleIndex= &pi;
//     filterArgument.rangeList = &rl;
//     filterArgument.siteDepthMin = (int)argSiteDepthMin.AsInteger();
//     filterArgument.siteDepthMax = (int)argSiteDepthMax.AsInteger();
//     filterArgument.indvDepthMin = (int)argIndvDepthMin.AsInteger();
//     filterArgument.indvDepthMax = (int)argIndvDepthMax.AsInteger();
//     filterArgument.indvQualMin = (int)argIndvQualMin.AsInteger();
//     filterArgument.infoGrepper = &ig;

//     VCFFile vcfFile;
//     vcfFile.open(argInputVCFFileName.c_str());
//     vcfFile.output012File(argOutputPrefix.c_str(), filterArgument);
//     vcfFile.close();

//     return 0;
// }

