#ifndef _VCFINPUTFILE_H_
#define _VCFINPUTFILE_H_

#include "tabix.h"
#include "VCFRecord.h"

#include "Logger.h"
#include "VCFFilter.h"

class VCFInputFile{
public:
    VCFInputFile (const char* fn):
        fp(NULL), tabixHandle(NULL){
        this->fileName = fn;
        // open file
        this->fp = new LineReader(fn);
        
        // read header
        while (this->fp->readLine(&_line)){
            if (_line[0] == '#') {
                this->header.push_back(_line);
                if (_line.substr(0, 6) == "#CHROM") {
                    this->record.createIndividual(_line);
                    this->headerLoaded = true;
                    break;
                }
                continue;
            }
            if (_line[0] != '#') {
                FATAL("Wrong VCF header");
            }
        }

        this->tabixHandle = 0;
        this->iter = 0;
        this->hasIndex = this->openIndex();

        this->rangeIdx = 0;
        this->clearRange();
        this->readOnlyLine = 0;

        this->lineLength = 1024;
        this->line = new char[1024];
        assert(this->line);

        /* this->mutableLine = (char*) malloc(sizeof(char) * 1024); */
        /* this->mutableLineLength = 1024; */
    };
    ~VCFInputFile(){
        closeIndex();
        this->record.deleteIndividual();
        if (this->fp) delete this->fp;
        /* if (this->mutableLine) { */
        /*     free(this->mutableLine); */
        /*     this->mutableLine = NULL; */
        /* } */
    };

    // use current subset of included people
    // to reconstruct a new VCF header line
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
			if (ti_lazy_index_load(this->tabixHandle) != 0) {
                // failed to open tabix index
                this->hasIndex = false;
                return false;
            } else {
                this->hasIndex = true;
                return true;
            }
        }
        return true;
    };
    void closeIndex(){
        if (this->iter) {
            ti_iter_destroy(this->iter);
            this->iter = 0;
        }
        ti_close(this->tabixHandle);
        this->tabixHandle = 0;
    };

    void clearRange() {
        this->range.clear();
        this->rangeIdx = 0;
        this->readOnlyLine = 0;
    };
    void setRangeFile(const char* fn) {
        if (!fn || strlen(fn) == 0)
            return;
        if (this->range.size()) {
            fprintf(stdout, "Clear existing %d regions.\n", this->range.size());
            this->clearRange();
        }
        this->range.addRangeFile(fn);
    };
    // @param l is a string of range(s)
    void setRange(const char* chrom, int begin, int end) {
        if (this->range.size()) {
            fprintf(stdout, "Clear existing %d regions.\n", this->range.size());
            this->clearRange();
        }
        this->range.addRange(chrom, begin, end);
    };
    void setRangeList(const char* l){
        if (!l || strlen(l) == 0)
            return;

        if (this->range.size()) {
            fprintf(stdout, "Clear existing %d regions.\n", this->range.size());
            this->clearRange();
        }
        this->range.addRangeList(l);
    };
    void setRangeList(RangeList& rl){
        if (rl.size() == 0) return;

        if (this->range.size()) {
            fprintf(stdout, "Clear existing %d regions.\n", this->range.size());
            this->clearRange();
        }
        this->range = rl;
    };
    /**
     * @return true: a valid VCFRecord 
     */
    bool readRecord(){
        assert(this->headerLoaded);
        
        // load contents 
        if (this->range.size() > 0) {
            if (this->hasIndex) {                 // there is index
                while (rangeIdx < this->range.size()) {
                    if (!this->readOnlyLine) { // last time does not read a valid line
                        // get range
                        this->range.obtainRange(rangeIdx, &_r);
                        // parse range
                        if (ti_parse_region(tabixHandle->idx, _r.c_str(), &_tid, &_beg, &_end) != 0){
                            FATAL("Cannot ti_parse_region");
                        }
                        iter =  ti_queryi(tabixHandle, _tid, _beg, _end);
                        this->readOnlyLine = ti_read(this->tabixHandle, iter, &_len);
                        if (this->readOnlyLine) { // s is valid
                            duplicateReadOnlyLine(this->readOnlyLine, _len);
                            this->record.parse(this->line, this->lineLength);
                            return true;
                        } else{ 
                            // continue to next rangeIdx
                            ti_iter_destroy(this->iter);
                            this->iter = 0;
                            rangeIdx ++;
                            continue;
                        }
                    } else {  // last time read a valid line
                        this->readOnlyLine = ti_read(this->tabixHandle, iter, &_len);
                        if (!this->readOnlyLine) {
                            rangeIdx ++;
                            continue;
                        } else {
                            duplicateReadOnlyLine(this->readOnlyLine, _len);
                            this->record.parse(this->line, this->lineLength);
                            return true;
                        }
                    }
                } // end while 
                return false;
            } // end hasIndex
            else { // no index
                // no index so force quitting
                fprintf(stderr, "Failed to load index when accessing VCF by region\n");
                abort();

                /* // legacy code (when there is no index but still want to parse VCF file). */
                /* // legacy code that can handle region in a slow but working way. */
                /* while(this->fp->readLine(&this->line)){ */
                /*     this->record.parse(line.c_str()); */
                /*     if (!this->range.isInRange(this->record.getChrom(), this->record.getPos())) */
                /*         continue; */
                /* } */
                /* this->record.parse(this->line.c_str()); */
            };
        } else { // go line by line
            if (this->fp->readLine(&_line)){
                duplicateReadOnlyLine(_line);
                this->record.parse(this->line, this->lineLength);
                return true;
            } else {
                return false; // reach file end
            }
        }
    };
    void includePeople(const char* s) {
        this->record.includePeople(s);
    };
    void includePeople(const std::vector<std::string>& v){
        this->record.includePeople(v);
    };
    void includePeopleFromFile(const char* fn) {
        this->record.includePeopleFromFile(fn);
    };
    void includeAllPeople(){
        this->record.includeAllPeople();
    };
    void excludePeople(const char* s) {
        this->record.excludePeople(s);
    };
    void excludePeople(const std::vector<std::string>& v){
        this->record.excludePeople(v);
    };
    void excludePeopleFromFile(const char* fn) {
        this->record.excludePeopleFromFile(fn);
    };
    void excludeAllPeople(){
        this->record.excludeAllPeople();
    };
    VCFRecord& getVCFRecord() {return this->record;};
    const char* getLine() const {return this->line;};
  private:
    // disable copy-constructor
    VCFInputFile(const VCFInputFile& v){};
    VCFInputFile& operator=(const VCFInputFile& v){};
    
    // copy this->readOnlyLine to this->mutableLine
    void duplicateReadOnlyLine(const char* s, int len){
        while (this->lineLength < len) {
            this->lineLength *= 2;
        }
        delete[] this->line;
        this->line = new char[this->lineLength];
        assert(this->line);
        strncpy(this->line, s, len + 1);
    };
    void duplicateReadOnlyLine(const std::string& s){
        int len = s.size();
        while (this->lineLength < len) {
            this->lineLength *= 2;
        }
        delete[] this->line;
        this->line = new char[this->lineLength];
        assert(this->line);
        strncpy(this->line, s.c_str(), len + 1);
    };
  private:
    VCFHeader header;
    VCFRecord record;
    VCFFilter filter;
    
    LineReader* fp;
    tabix_t * tabixHandle;
    RangeList range;
    
    std::string fileName;

    bool headerLoaded;
    bool hasIndex;

    // variable used for accessing by region.
    int rangeIdx;
    ti_iter_t iter; 
    char* line;
    int lineLength;           // line length (excluding the trailing \n)
    // std::string line;         // actual data line (only used for iterative parsing VCF).
    const char* readOnlyLine; // read-only copy from tabix
    
    // those variables are declared for speed up readRecord()
    std::string _r;       // region
    int _tid, _beg, _end;   // rgeion for tabix
    int _len; 
    std::string _line;

    /* char* mutableLine;        // will contain cached copy of readOnlyLine to sparse up parsing. */
    /* int mutableLineLength;     */

};

#endif /* _VCFINPUTFILE_H_ */
