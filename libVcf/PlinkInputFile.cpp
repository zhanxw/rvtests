#include "PlinkInputFile.h"
#include "SimpleMatrix.h"

// @param m: people by marker matrix
int PlinkInputFile::readIntoMatrix(SimpleMatrix* mat) {
    assert(mat);

    // read bed
    int numPeople = this->getNumIndv();
    int numMarker = this->getNumMarker();
    fprintf(stderr, "%d people, %d marker\n", numPeople, numMarker);
    if (snpMajorMode) {
        // unsigned char mask[] = { 0x3, 0xc, 0x30, 0xc0 }; //0b11, 0b1100, 0b110000, 0b11000000
        unsigned char mask = 0x3; // 0b0000011
        unsigned char c;
        (*mat).resize( numPeople, numMarker );
        for (int m = 0; m < numMarker; m++){
            for (int p = 0; p < numPeople; p++) {
                int offset = p & (4 - 1);
                if (offset == 0) {
                    int ret = fread(&c, sizeof(unsigned char), 1, fpBed);
                    assert (ret == 1);
                }
                unsigned char geno = (c  >> (offset << 1)) & mask;
                switch (geno){
                    case HOM_REF:
                        (*mat)[p][m] = 0;
                        break;
                    case HET:
                        (*mat)[p][m] = 1;
                        break;
                    case HOM_ALT:
                        (*mat)[p][m] = 2;
                        break;
                    case MISSING:
                        (*mat)[p][m] = -9;
                        break;
                    default:
                        REPORT("Read PLINK genotype error!\n");
                        break;
                };
                //
                //unsigned char temp =  (c >> (offset << 1) ) & 3 ;
                //printf("m=%d, p=%d, offset=%d, c=%d, geno=%d, geno = %d, temp=%d\n", m,p,offset,c, geno, (int)(*mat)[p][m],temp);
            }
        }
    } else { // Indv_Major_Mode
        unsigned char mask[] = { 0x3, 0xc, 0x30, 0xc0 }; //0b11, 0b1100, 0b110000, 0b11000000
        unsigned char c;
        (*mat).resize( numPeople, numMarker );
        for (int p = 0; p < numPeople; p++) {
            for (int m = 0; m < numMarker; m++){
                int offset = m & (4 - 1);
                if (offset == 0) {
                    int ret = fread(&c, sizeof(unsigned char), 1, fpBed);
                    assert( ret == 1);
                }
                unsigned char geno = (c & mask[offset]) >> (offset << 1);
                switch (geno){
                    case HOM_REF:
                        (*mat)[m][p] = 0;
                        break;
                    case HET:
                        (*mat)[m][p] = 1;
                        break;
                    case HOM_ALT:
                        (*mat)[m][p] = 2;
                        break;
                    case MISSING:
                        (*mat)[m][p] = -9;
                        break;
                    default:
                        REPORT("Read PLINK genotype error!\n");
                        break;
                };
            }
        }
    }
    return this->getNumMarker() * this->getNumIndv();
};

int PlinkInputFile::readIntoMatrix(SimpleMatrix* mat, std::vector<std::string>* peopleNames, std::vector<std::string>* markerNames) {
    assert (mat);

    // read bed
    int numPeople = getNumIndv();
    int numMarker = getNumMarker();

    mat->resize( peopleNames == NULL ? numPeople: peopleNames->size(),
                    markerNames == NULL ? numMarker: markerNames->size());

    // get people index
    std::vector<int> peopleIdx;
    if (peopleNames == NULL || peopleNames->size() == 0) {
        // all peoples
        peopleIdx.resize(pid2Idx.size());
        for (unsigned int i = 0; i < pid2Idx.size(); i++)
            peopleIdx[i] = (i);
    } else {
        for (unsigned int i = 0; i < peopleNames->size(); i++)
            if ( pid2Idx.find( (*peopleNames)[i] ) != pid2Idx.end()) {
                peopleIdx.push_back(  pid2Idx[ (*peopleNames)[i] ]) ;
            }
    }

    // get marker index
    std::vector<int> markerIdx;
    if (markerNames == NULL || markerNames->size() == 0) {
        // all markers
        markerIdx.resize(snp.size());
        for (unsigned int i = 0; i < snp.size(); i++)
            markerIdx[i]= (i);
    } else {
        for (unsigned int i = 0; i < markerNames->size(); i++)
            if ( snp2Idx.find( (*markerNames)[i] ) != snp2Idx.end()) {
                markerIdx.push_back(  snp2Idx[ (*markerNames)[i] ]) ;
            }
    }

    int peopleToRead = peopleIdx.size();
    int markerToRead = markerIdx.size();

    unsigned char mask[] = { 0x3, 0xc, 0x30, 0xc0 }; //0b11, 0b1100, 0b110000, 0b11000000
    (*mat).resize( peopleToRead, markerToRead);
    if (snpMajorMode) {
        for (int p = 0; p < peopleToRead; p++) {
            for (int m = 0; m < markerToRead; m++) {
                // get file position
                int pos = 3 + (numPeople / 4 + 1 ) * markerIdx[m] + peopleIdx[p] / 4;
                int offset = peopleIdx[p] % 4;
                unsigned char c;
                fseek(this->fpBed, pos, SEEK_SET);
                fread(&c, sizeof(unsigned char), 1, fpBed);
                unsigned char geno = (c & mask[offset]) >> (offset << 1);
                switch (geno){
                    case HOM_REF:
                        (*mat)[p][m] = 0;
                        break;
                    case HET:
                        (*mat)[p][m] = 1;
                        break;
                    case HOM_ALT:
                        (*mat)[p][m] = 2;
                        break;
                    case MISSING:
                        (*mat)[p][m] = -9;
                        break;
                    default:
                        REPORT("Read PLINK genotype error!\n");
                        break;
                };
            }
        }
    } else { // Indv_Major_Mode
        for (int p = 0; p < numPeople; p++) {
            for (int m = 0; m < numMarker; m++){
                // get file position
                int pos = 3 + (numMarker / 4  + 1) * peopleIdx[p] + markerIdx[m] / 4;
                int offset = markerIdx[m] % 4;
                unsigned char c;
                fseek(this->fpBed, pos, SEEK_SET);
                fread(&c, sizeof(unsigned char), 1, fpBed);

                unsigned char geno = (c & mask[offset]) >> (offset << 1);
                switch (geno){
                    case HOM_REF:
                        (*mat)[m][p] = 0;
                        break;
                    case HET:
                        (*mat)[m][p] = 1;
                        break;
                    case HOM_ALT:
                        (*mat)[m][p] = 2;
                        break;
                    case MISSING:
                        (*mat)[m][p] = -9;
                        break;
                    default:
                        REPORT("Read PLINK genotype error!\n");
                        break;
                };
            }
        }
    }
    return this->getNumMarker() * this->getNumIndv();    
};
