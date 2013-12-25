#include "VCFUtil.h"

int main() {
  {
    // read all sites
    VCFInputFile vin("test.vcf.gz");
    int lineNo = 0;
    while (vin.readRecord()){
        lineNo ++;
        VCFRecord& r = vin.getVCFRecord(); 
        VCFPeople& people = r.getPeople();
        VCFIndividual* indv;

        printf("%s:%d\n", r.getChrom(), r.getPos());
    };
    fprintf(stdout, "Total %d VCF records have converted successfully\n", lineNo);
  }
  
  {
    // read sites
    VCFInputFile vin("test.vcf.gz");
    vin.setSiteFile("testVCFSite.sites");
    int lineNo = 0;
    while (vin.readRecord()){
        lineNo ++;
        VCFRecord& r = vin.getVCFRecord(); 
        VCFPeople& people = r.getPeople();
        VCFIndividual* indv;

        printf("%s:%d\n", r.getChrom(), r.getPos());
    };
    fprintf(stdout, "Total %d VCF records have converted successfully\n", lineNo);
  }

  {
    // read sites within given range
    VCFInputFile vin("test.vcf.gz");
    vin.setRangeList("1:196341181-196341799");
    vin.setSiteFile("testVCFSite.sites");    
    int lineNo = 0;
    while (vin.readRecord()){
        lineNo ++;
        VCFRecord& r = vin.getVCFRecord(); 
        VCFPeople& people = r.getPeople();
        VCFIndividual* indv;

        printf("%s:%d\n", r.getChrom(), r.getPos());
    };
    fprintf(stdout, "Total %d VCF records have converted successfully\n", lineNo);
  }
};
