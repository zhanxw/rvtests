#include "VCFInputFile.h"
#include "Utils.h"
#include "IO.h"
// use current subset of included people
// to reconstruct a new VCF header line
void VCFInputFile::rewriteVCFHeader() {
  std::string s = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  VCFPeople& people = this->record.getPeople();
  for (unsigned int i = 0; i <people.size(); i++ ){
    s += '\t';
    s += people[i]->getName();
  }
  this->header[this->header.size()-1] = s;
};

void VCFInputFile::setRangeMode() {
  if (!this->hasIndex){
    fprintf(stderr, "Please provide VCF index when accessing VCF by region\n");
    abort();
  };
  this->mode = VCFInputFile::RANGE_MODE;
  this->rangeBegin = this->range.begin();
  this->rangeEnd = this->range.end();
  this->rangeIterator = this->range.begin();
  
};
void VCFInputFile::clearRange() {
#ifndef NDEBUG
  if (this->range.size()) {
    fprintf(stderr, "Clear existing %zu range.\n", this->range.size());
  }
#endif
  this->range.clear();
  this->ti_line = 0;
};

/**
 * @param fn: the file contains two column: old_id new_id
 */
int VCFInputFile::updateId(const char* fn){
  // load conversion table
  LineReader lr(fn);
  std::map<std::string, std::string> tbl;
  std::vector<std::string> fd;
  while(lr.readLineBySep(&fd, "\t ")){
    if (tbl.find(fd[0]) != tbl.end()) {
      fprintf(stderr, "Duplicated original ids: [ %s ], replace it to new id anyway.\n", fd[0].c_str());
    };
    tbl[fd[0]] = fd[1];
  }

  // rewrite each people's name
  std::string s = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  VCFPeople& people = this->record.getPeople();
  int n = 0;
  for (unsigned int i = 0; i <people.size(); i++ ){
    if (tbl.find(people[i]->getName()) != tbl.end()) {
      ++n;
      people[i]->setName(tbl[people[i]->getName()]);
    }
  }
  this->rewriteVCFHeader();

  // return result
  return n;

}
