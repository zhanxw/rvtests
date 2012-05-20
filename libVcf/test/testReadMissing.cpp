#include "VCFUtil.h"

int main() {
  {
    puts("Check whether field is missing for each individual");
    VCFInputFile vin("missing.vcf.gz");
    while (vin.readRecord()){
      VCFRecord& r = vin.getVCFRecord();
      VCFPeople& people = r.getPeople();
      VCFIndividual* indv;

      bool missing; // missing indicator
      for (int i = 0; i < people.size(); i++){
        indv = people[i];
        for (int idx = 0; idx < 5; idx++) {
          const VCFValue& v = indv->get(idx, &missing);
          if (!missing) {
            printf(".");
          } else {
            printf("M");
          }
        }
        printf("\t");
      };

      printf("\n");
    };
    puts(". - Individual's column is not missing");
    puts("M - Individual's column is missing\n");

  }
  {
    puts("Check whether genotype is missing for each individual");
    VCFInputFile vin("missing.vcf.gz");
    while (vin.readRecord()){
      VCFRecord& r = vin.getVCFRecord();
      VCFPeople& people = r.getPeople();
      VCFIndividual* indv;

      bool missingGenotype; // missing indicator
      int GTidx = r.getFormatIndex("GT");
      for (int i = 0; i < people.size(); i++){
        indv = people[i];
        const VCFValue& gt = indv->justGet(GTidx);
        missingGenotype = gt.isMissingGenotype();
        if (!missingGenotype) {
          printf(".");
        } else {
          printf("M");
        }
        printf("\t");
      }
      printf("\n");
    }
    puts(". - Individual's genotype column is not missing");
    puts("M - Individual's genotype column is missing\n");
  }
};
