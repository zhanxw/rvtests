#ifndef _VCFRECORD_H_
#define _VCFRECORD_H_

#include "VCFFunction.h"
#include "VCFIndividual.h"
#include "VCFInfo.h"
#include "OrderedMap.h"

typedef OrderedMap<int, VCFIndividual*> VCFPeople;

class VCFRecord{
public:
    VCFRecord(){
        this->hasAccess = false;
    };

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
            FATAL("not enough people in the VCF (VCF does not contain genotype and individuals?)");
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

    //////////////////////////////////////////////////////////////////////
    // Code related with include/exclude people
    void includePeople(const std::string& name){
        if (name.size() == 0) return;
        for (unsigned int i = 0 ; i != this->allIndv.size() ; i++) {
            VCFIndividual* p = this->allIndv[i];
            if (p->getName() == name) {
                p->include();
            }
        }
        this->hasAccess = false;
    };
    void includePeople(const std::vector<std::string>& v){
        for (unsigned int i = 0; i < v.size(); i++){
            this->includePeople(v[i]);
        }
        this->hasAccess = false;
    };
    void includePeopleFromFile(const char* fn){
        if (!fn || strlen(fn) == 0) return;
        LineReader lr(fn);
        std::vector<std::string> fd;
        while(lr.readLineBySep(&fd, "\t ")) {
            for (unsigned int i = 0; i < fd.size(); i++)
                this->includePeople(fd[i]);
        }
        this->hasAccess = false;
    };
    void includeAllPeople() {
        for (unsigned int i = 0 ; i != this->allIndv.size() ; i++) {
            VCFIndividual* p = this->allIndv[i];
            p->include();
        }
        this->hasAccess = false;
    };
    void excludePeople(const std::string& name){
        if (name.size() == 0) return;
        for (unsigned int i = 0 ; i != this->allIndv.size() ; i++) {
            VCFIndividual* p = this->allIndv[i];
            if (p->getName() == name) {
                p->exclude();
            }
        }
        this->hasAccess = false;
    };
    void excludePeople(const std::vector<std::string>& v){
        for (unsigned int i = 0; i != v.size(); i++ ){
            this->excludePeople(v[i]);
        }
        this->hasAccess = false;
    };
    void excludePeopleFromFile(const char* fn){
        if (!fn || strlen(fn) == 0) return;
        LineReader lr(fn);
        std::vector<std::string> fd;
        while(lr.readLineBySep(&fd, "\t ")) {
            for (unsigned int i = 0; i != fd.size(); i++)
                this->excludePeople(fd[i]);
        }
        this->hasAccess = false;
    };
    void excludeAllPeople() {
        for (unsigned int i = 0 ; i != this->allIndv.size() ; i++) {
            VCFIndividual* p = this->allIndv[i];
            p->exclude();
        }
        this->hasAccess = false;
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
    // const char* getLine() {return this->line;};
    VCFPeople& getPeople(){
        if (!this->hasAccess) {
            this->selectedIndv.clear();
            for (int i = 0; i < this->allIndv.size(); i++){
                if (allIndv[i]->isInUse()) {
                    this->selectedIndv[this->selectedIndv.size()] = allIndv[i];
                }
            }
            this->hasAccess = true;
        }
        return this->selectedIndv;
    };
    int getFormatIndex(const char* s){
        //TODO: cache query
        int b = this->format.beg;
        int e = this->format.end;
        int l = strlen(s);
        int idx = 0;

        // locate first field
        while ( b < e) {
            // check match
            bool match = true;
            for (int i = 0; i < l ; i++) {
                if (line[b + i]  != s[i]) {
                    match = false;
                    break;
                }
            }
            if (match) return idx;
            else {
                // skip to next field
                idx++;
                while ( line[b++] != ':') {
                    if (b >= e) {
                        return -1;
                    }
                }
            }
        }
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

    // indicates if getPeople() has been called
    bool hasAccess;

}; // VCFRecord

#endif /* _VCFRECORD_H_ */
