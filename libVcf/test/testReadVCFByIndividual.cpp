#include "VCFUtil.h"

int main() {
  VCFInputFile vin("test.vcf");
  vin.excludeAllPeople();
  vin.includePeople("111409786");
  vin.includePeople("113539618");
  int lineNo = 0;
  while (vin.readRecord()){
    lineNo ++;
    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    printf("%s:%d\t", r.getChrom(), r.getPos());

    // e.g.: get TAG from INFO field
    // fprintf(stderr, "%s\n", r.getInfoTag("ANNO"));

    // e.g.: Loop each (selected) people in the same order as in the VCF
    for (int i = 0; i < people.size(); i++) {
      indv = people[i];
      // get GT index. if you are sure the index will not change, call this function only once!
      int GTidx = r.getFormatIndex("GT");
      if (GTidx >= 0)
        printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first field of each individual
      else
        fprintf(stderr, "Cannot find GT field!\n");
    }
    printf("\n");
  };
  fprintf(stdout, "Total %d VCF records have converted successfully\n", lineNo);

};
