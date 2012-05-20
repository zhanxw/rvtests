#include "VCFUtil.h"

int main() {
    VCFInputFile vin("missing.vcf.gz");
    vin.setRangeList("1:0");
    while (vin.readRecord()){
        VCFRecord& r = vin.getVCFRecord(); 
        printf("%s:%d\t", r.getChrom(), r.getPos());
        printf("\n");
    };

    vin.setRangeList("X:800-900");
    while (vin.readRecord()){
        VCFRecord& r = vin.getVCFRecord(); 
        printf("%s:%d\t", r.getChrom(), r.getPos());
        printf("\n");
    };

    vin.setRangeList("1:30-50");
    while (vin.readRecord()){
        VCFRecord& r = vin.getVCFRecord(); 
        printf("%s:%d\t", r.getChrom(), r.getPos());
        printf("\n");
    };
};
