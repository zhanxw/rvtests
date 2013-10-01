#include "VCFUtil.h"

int main() {
  VCFInputFile vin("test.bcf.gz");
  vin.setRangeList("1:196341364-196341449");
  while (vin.readRecord()){
    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;
    
    printf("%s:%d\t", r.getChrom(), r.getPos());

    bool tagMissing;
    VCFInfo& info = r.getVCFInfo();
    std::string anno = info.getTag("ANNO", &tagMissing).toStr();
    
    printf("ANNO=%s\t", anno.c_str());
    // assert(tagMissing); // all variant has tagMissing == true

    bool missingGenotype; // missing indicator
    int GTidx = r.getFormatIndex("GT");
    for (int i = 0; i < people.size(); i++){
      indv = people[i];
      const VCFValue& gt = indv->justGet(GTidx);
      missingGenotype = gt.isMissingGenotype();
      if (!missingGenotype) {
        printf("%d", gt.getGenotype());
      } else {
        printf("M");
      }
      printf("\t");
    }
    printf("\n");
  }

  return 0;
}
