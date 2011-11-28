#ifndef _VCFINPUTFILE_H_
#define _VCFINPUTFILE_H_

#include "tabix.h"
#include "VCFRecord.h"

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

        this->rangeIdx = 0;
        this->s = 0;
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
    // @param l is a string of range(s)
    void setRangeList(const char* l){
        this->range.addRangeList(l);
    }

    bool readRecord(){
        assert(this->headerLoaded);
        // load contents 
        if (this->range.size() > 0) {
            if (this->hasIndex) {                 // there is index
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

    // variable used for accessing by region.
    int rangeIdx;
    ti_iter_t iter; 
    const char* s;

};

#endif /* _VCFINPUTFILE_H_ */
