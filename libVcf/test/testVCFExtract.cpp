#include "VCFUtil.h"

void print(VCFInputFile* v) {
  VCFInputFile& vin = *v;
    int lineNo = 0;
    printf("CHROM:POS\tAC\tAN\tAF\tGENOTYPES...\n");
    while (vin.readRecord()){
        lineNo ++;
        VCFRecord& r = vin.getVCFRecord(); 
        VCFPeople& people = r.getPeople();
        VCFIndividual* indv;

        bool missing;
        printf("%s:%d\t%s\t%s\t%s\t",
               r.getChrom(), r.getPos(),
               r.getInfoTag("AC", &missing).toStr(),
               r.getInfoTag("AN", &missing).toStr(),
               r.getInfoTag("AF", &missing).toStr()
               );
        
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
    
};
int main(){
  {
    VCFInputFile vin("test.vcf");
    print(&vin);
    vin.close();
  }
  puts("--------- calculated AF > 0.001 -------------------------------------------------------------");
  {
    VCFExtractor ve("test.vcf");
    ve.setSiteFreqMin(0.001);
    print(&ve);
    ve.close();
  }
  puts("--------- AF > 0.001 -------------------------------------------------------------");
  {
    VCFExtractor ve("test.vcf");
    ve.setSiteFreqMin(0.001);
    ve.setUseSiteFreqFromInfo();
    print(&ve);
    ve.close();
  }
  puts("--------- Variant Site Only -------------------------------------------------------------");
  {
    VCFExtractor ve("test.vcf");
    ve.setVariantSiteOnly();
    print(&ve);
    ve.close();
  }
  
  return 0;
};
