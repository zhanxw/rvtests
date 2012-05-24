#include "VCFUtil.h"

int main() {
  VCFInputFile vin("test.vcf");
  int lineNo = 0;
  while (vin.readRecord()){
    lineNo ++;
    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;


    VCFInfo& info = r.getVCFInfo();

    // output HM3, DP, CBR, ANNO and FAKE

    // (1) loop to access tag
    bool inHapMap = false;    
    int dp;
    float cbr;
    std::string anno;
    std::string fake;
    
    std::string tagName;
    VCFValue value;
    for (unsigned int i = 0; i < info.size() ; i++) {
      info.at(i, &tagName, &value);
      if (tagName == "HM3") {
        inHapMap = true;
      } else if (tagName == "DP") {
        dp = value.toInt();
      } else if (tagName == "CBR") {
        cbr = value.toDouble();
      } else if (tagName == "ANNO") {
        anno = value.toStr();
      } else if (tagName == "FAKE") {
        fake = value.toStr();
      }
    }
    assert(fake == "");
    printf( "%s:%d\t", r.getChrom(), r.getPos());
    printf( "%s\t%d\t%f\t%s\n", inHapMap?"inHapmap":"notHapMap", dp, cbr, anno.c_str());
        
    // (2) access by tag
    bool tagMissing;
    info.getTag("HM3", &tagMissing);
    inHapMap = tagMissing;
    dp       = info.getTag("DP", &tagMissing).toInt();
    cbr      = info.getTag("CBR", &tagMissing).toDouble();
    anno     = info.getTag("ANNO", &tagMissing).toStr();
    fake     = info.getTag("FAKE", &tagMissing).toStr();
    printf( "%s:%d\t", r.getChrom(), r.getPos());
    assert(fake == "" && tagMissing == true);

    printf( "%s\t%d\t%f\t%s\n", inHapMap?"inHapmap":"notHapMap", dp, cbr, anno.c_str());

    // output results


    printf("\n");
  };
  printf( "Total %d VCF records have converted successfully\n", lineNo);

};
