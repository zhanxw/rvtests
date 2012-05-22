#include "VCFInputFile.h"

// use current subset of included people
// to reconstruct a new VCF header line
void VCFInputFile::rewriteVCFHeader() {
  std::string s = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  VCFPeople& people = this->record.getPeople();
  for (int i = 0; i <people.size(); i++ ){
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
  if (this->range.size()) {
    fprintf(stdout, "Clear existing %zu range.\n", this->range.size());
  }
  this->range.clear();
  this->ti_line = 0;
};
