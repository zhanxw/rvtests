#ifndef _VCFRECORD_H_
#define _VCFRECORD_H_

#include "VCFFunction.h"
#include "VCFRecord.h"
#include "VCFIndividual.h"
#include "VCFInfo.h"
#include "OrderedMap.h"

#if 0 
extern int parseTillChar(const char c, const char* line, const int beg, VCFValue* ret); 
extern int parseTillChar(const char* c, const char* line, const int beg, VCFValue* ret);
#endif 

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

#endif /* _VCFRECORD_H_ */
