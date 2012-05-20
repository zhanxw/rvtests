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

  void parse(char* line, int len){
    this->vcfInfo.reset();
    this->self.line = line;
    this->self.beg = 0;
    this->self.end = len;

    // this->line = line;
    // go through VCF sites (first 9 columns)

    assert(this->self.parseTill('\t', &this->chrom) == 0);
    line[this->chrom.end] = '\0';

    assert(this->self.parseTill('\t', this->chrom.end + 1, &this->pos) == 0);
    line[this->pos.end] = '\0';

    assert(this->self.parseTill('\t', this->pos.end + 1, &this->id) == 0);
    line[this->id.end] = '\0';

    assert(this->self.parseTill('\t', this->id.end + 1, &this->ref) == 0);
    line[this->ref.end] = '\0';

    assert(this->self.parseTill('\t', this->ref.end + 1, &this->alt) == 0);
    line[this->alt.end] = '\0';

    assert(this->self.parseTill('\t', this->alt.end + 1, &this->qual) == 0);
    line[this->qual.end] = '\0';

    assert(this->self.parseTill('\t', this->qual.end + 1, &this->filt) == 0);
    line[this->filt.end] = '\0';

    assert(this->self.parseTill('\t', this->filt.end + 1, &this->info) == 0);
    line[this->info.end] = '\0';
    this->vcfInfo.parse(this->info); // lazy parse inside VCFInfo

    assert(this->self.parseTill('\t', this->info.end + 1, &this->format) == 0);
    line[this->format.end] = '\0';

    // now comes each individual genotype
    int idx = 0; // peopleIdx
    VCFIndividual* p = this->allIndv[idx];
    int end = this->format.end;
    VCFValue indv;
    while( this->self.parseTill('\t', end + 1, &indv) == 0) {
      p->parse(indv);
      end = indv.end;

      if (indv.end == this->self.end) {
        fprintf(stderr, "indv.end = %d, self.end = %d\n", indv.end, this->self.end);
        break;
      }
      // check next individual
      idx ++ ;
      if (idx >= this->allIndv.size()){
        fprintf(stderr, "Idx = %d, allIndv.size() = %d \n", idx, this->allIndv.size());
        FATAL("VCF header have LESS people than VCF content!");
      }
      p = this->allIndv[idx];
    };

    if (idx + 1 > this->allIndv.size()) {
      FATAL("VCF header have MORE people than VCF content!");
    };
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
  VCFInfo& getVCFInfo() {
    return this->vcfInfo;
  };
  const VCFValue& getInfoTag(const char* tag, bool* exists) {
    return this->vcfInfo.getTag(tag, exists);
  };
  void output(FILE* fp) const{
    this->self.output(fp);
    fputc('\n', fp);
  };
public:
  const char* getChrom() { return this->chrom.toStr(); };
  const int   getPos()  { return this->pos.toInt(); };
  const char* getID() { return this->id.toStr(); };
  const char* getRef() { return this->ref.toStr(); };
  const char* getAlt() { return this->alt.toStr(); };
  const char* getQual() { return this->qual.toStr(); };
  const char* getFilt() { return this->filt.toStr(); };
  const char* getInfo() { return this->vcfInfo.toStr(); };
  const char* getFormat() { return this->format.toStr(); };

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
  /**
   * You want to know tag "GQ" where FORMAT column is "GT:GQ:PL"
   * then call getFormatIndx("GQ") will @return 1 (0-based index)
   * @return -1 when not found
   */
  int getFormatIndex(const char* s){
    int b = this->format.beg;
    int e = this->format.end;
    int idx = 0;

    // locate first field
    while ( b < e) {
      // check match
      bool match = true;
      for (int i = 0; s[i] != '\0' ; i++) {
        if (this->format.line[b + i]  != s[i]) {
          match = false;
          break;
        }
      }
      if (match) return idx;

      // skip to next field
      idx++;
      while ( this->format.line[b++] != ':') {
        if (b >= e) {
          return -1;
        }
      }
    }
  };
  const VCFValue& getSelf() const{
    return this->self;
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

  VCFInfo vcfInfo;

  VCFValue self;       // a self value points to itself
  // const char* line; // points to data line

  // indicates if getPeople() has been called
  bool hasAccess;

}; // VCFRecord

#endif /* _VCFRECORD_H_ */
