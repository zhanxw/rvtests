#include "PlinkOutputFile.h"
#include "SimpleMatrix.h"

void PlinkOutputFile::writeBED(SimpleMatrix* mat, int nPeople, int nMarker){
    /* int nPeople = mat->cols; */
    /* int nMarker = mat->rows; */
    unsigned char c = 0;
    int offset;
    for (int m = 0; m < nMarker; m++){
        for (int i = 0; i < nPeople ; i ++) {
            offset = i & (4 - 1);
            int geno = (int)( (*mat)[m][i]);
            switch(geno){
                case HOM_REF:
                    setGenotype(&c, offset, HET); // het: 0b01
                    break;
                case HET:
                    setGenotype(&c, offset, HET); // het: 0b01
                    break;
                case HOM_ALT:
                    setGenotype(&c, offset, HOM_ALT); // hom alt: 0b11
                    break;
                default:
                    setGenotype(&c, offset, MISSING); // missing
                    break;
            }
        }
        if ( offset == 3) { // 3: 4 - 1, so every 4 genotype we will flush
            fwrite(&c, sizeof(char), 1, this->fpBed);
            c = 0;
        }
    };
    if (nPeople % 4 != 0 )
        fwrite(&c, sizeof(char), 1, this->fpBed);
};
