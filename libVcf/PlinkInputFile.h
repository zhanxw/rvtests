#ifndef _PLINKINPUTFILE_H_
#define _PLINKINPUTFILE_H_

#include <map>
#include <string>

#include "Exception.h"
#include "IO.h"
#include "MathMatrix.h"


class PlinkInputFile{
public:
    PlinkInputFile(const char* fnPrefix) {
        this->prefix = fnPrefix;
        this->fpBed = fopen( (prefix + ".bed").c_str(), "rb");
        this->fpBim = fopen( (prefix + ".bim").c_str(), "rt");
        this->fpFam = fopen( (prefix + ".fam").c_str(), "rt");
        if (!this->fpBed || !this->fpBim || !this->fpFam){
            REPORT("Cannot open binary PLINK file!");
            abort();
        }
        // write Bed header
        char c;
        // magic number
        char magic1 = 0x6c; // 0b01101100;
        int ret;
        ret = fread(&c, sizeof(char), 1, this->fpBed);
        assert(ret == 1);
        if (c != magic1) {
            fprintf(stderr, "Magic number of binary PLINK file does not match!\n");
            abort();
        }
        int magic2 = 0x1b; // 0b00011011;
        ret = fread(&c, sizeof(char), 1, this->fpBed);
        assert(ret == 1);
        if (c != magic2) {
            fprintf(stderr, "Magic number of binary PLINK file does not match!\n");
            abort();
        }

        // snp major mode
        const int SNP_MAJOR_MODE = 0x01; //0b00000001;
        const int INDV_MAJOR_MODE = 0x00;
        ret = fread(&c, sizeof(char), 1, this->fpBed);
        assert(ret == 1);
        if ( c == SNP_MAJOR_MODE) {
            this->snpMajorMode = true;
        }else if ( c== INDV_MAJOR_MODE) {
            this->snpMajorMode = false;
        } else{
            fprintf(stderr, "Unrecognized major mode in binary PLINK file.\n");
            abort();
        }

        // read bim
        LineReader* lr = new LineReader( (this->prefix + ".bim").c_str());
        std::vector<std::string> fd;
        while( lr->readLineBySep(&fd, " \t") ){
            if (fd.size() != 6) {
                fprintf(stderr, "Wrong format in bim file.\n");
                continue;
            }

            if (snp2Idx.find(fd[1]) == snp2Idx.end()) {
                chrom.push_back (fd[0]);
                snp.push_back(fd[1]);
                snp2Idx[fd[1]] = snp2Idx.size() - 1;      //  -1 since fd[1] will be created using []
                mapDist.push_back(atof(fd[2].c_str()));
                pos.push_back(atoi(fd[3].c_str()));
                ref.push_back(fd[4][0]);
                alt.push_back(fd[5][0]);
            } else {
                fprintf(stderr, "duplicate marker name [%s], ignore!\n", fd[1].c_str());
            }
        };
        delete lr;

        // read fam
        lr = new LineReader( (this->prefix + ".fam").c_str());
        while( lr->readLineBySep(&fd, " \t") ){
            if (fd.size() != 6) {
                fprintf(stderr, "Wrong format in fam file.\n");
                continue;
            }

            // will skip loading fam, pid, mid
            if (pid2Idx.find(fd[1]) == pid2Idx.end()) {
                pid2Idx[fd[1]] = pid2Idx.size() - 1; //  -1 since fd[1] will be created using []
                indv.push_back(fd[1]);
                sex.push_back(atoi(fd[4].c_str()));
                pheno.push_back(atof(fd[5].c_str()));
            } else {
                fprintf(stderr, "duplicated person id [ %s ], ignore!\n", fd[1].c_str());
            }
        };
        delete lr;

    };
    ~PlinkInputFile() {
        fclose(this->fpBed);
        fclose(this->fpBim);
        fclose(this->fpFam);
    };

    // @param m: people by marker matrix
    int readIntoMatrix(Matrix* mat) {
        assert(mat);

        // read bed
        int numPeople = getNumIndv();
        int numMarker = getNumMarker();

        if (snpMajorMode) {
            // unsigned char mask[] = { 0x3, 0xc, 0x30, 0xc0 }; //0b11, 0b1100, 0b110000, 0b11000000
            unsigned char mask = 0x3; // 0b0000011
            unsigned char c;
            (*mat).Dimension( numPeople, numMarker );
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
            (*mat).Dimension( numPeople, numMarker );
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
    };


    int readIntoMatrix(Matrix* mat, std::vector<std::string>* peopleNames, std::vector<std::string>* markerNames) {
        assert (mat);

        // read bed
        int numPeople = getNumIndv();
        int numMarker = getNumMarker();

        mat->Dimension( peopleNames == NULL ? numPeople: peopleNames->size(),
                        markerNames == NULL ? numMarker: markerNames->size());

        // get people index
        std::vector<int> peopleIdx;
        if (peopleNames == NULL || peopleNames->size() == 0) {
            // all peoples
            peopleIdx.resize(pid2Idx.size());
            for (int i = 0; i < pid2Idx.size(); i++)
                peopleIdx[i] = (i);
        } else {
            for (int i = 0; i < peopleNames->size(); i++)
                if ( pid2Idx.find( (*peopleNames)[i] ) != pid2Idx.end()) {
                    peopleIdx.push_back(  pid2Idx[ (*peopleNames)[i] ]) ;
                }
        }

        // get marker index
        std::vector<int> markerIdx;
        if (markerNames == NULL || markerNames->size() == 0) {
            // all markers
            markerIdx.resize(snp.size());
            for (int i = 0; i < snp.size(); i++)
                markerIdx[i]= (i);
        } else {
            for (int i = 0; i < markerNames->size(); i++)
                if ( snp2Idx.find( (*markerNames)[i] ) != snp2Idx.end()) {
                    markerIdx.push_back(  snp2Idx[ (*markerNames)[i] ]) ;
                }
        }

        int peopleToRead = peopleIdx.size();
        int markerToRead = markerIdx.size();

        unsigned char mask[] = { 0x3, 0xc, 0x30, 0xc0 }; //0b11, 0b1100, 0b110000, 0b11000000
        (*mat).Dimension( peopleToRead, markerToRead);
        if (snpMajorMode) {
            for (int p = 0; p < peopleToRead; p++) {
                for (int m = 0; m < markerToRead; m++) {
                    // get file position
                    int pos = 3 + (numPeople / 4 + 1 ) * markerIdx[m] + peopleIdx[p] / 4;
                    int offset = peopleIdx[p] % 4;
                    unsigned char c;
                    fseek(this->fpBed, pos, SEEK_SET);
                    int ret = fread(&c, sizeof(unsigned char), 1, fpBed);                    
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
                    int ret = fread(&c, sizeof(unsigned char), 1, fpBed);

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
    };
    int getMarkerIdx(const std::string& m){
        if (this->snp2Idx.find(m) == this->snp2Idx.end()){
            return -1;
        } else{
            return (this->snp2Idx[m]);
        }
    };

    int getNumIndv() const {
        this->indv.size();
    };
    int getNumMarker() const {
        this->snp2Idx.size();
    };
public:
    std::vector<std::string> chrom;
    std::vector<std::string> snp;
    std::vector<double> mapDist;
    std::vector<int> pos;
    std::vector<char> ref;
    std::vector<char> alt;

    std::vector<std::string> indv; /// people ids
    std::vector<int> sex;
    std::vector<double> pheno;

private:
    std::map<std::string, int> snp2Idx;
    std::map<std::string, int> pid2Idx;
    // we reverse the two bits as defined in PLINK format,
    // so we can process 2-bit at a time.
    const static unsigned char HOM_REF = 0x0;     //0b00;
    const static unsigned char HET = 0x2;         //0b10;
    const static unsigned char HOM_ALT = 0x3;     //0b11;
    const static unsigned char MISSING = 0x1;     //0b01;

    FILE* fpBed;
    FILE* fpBim;
    FILE* fpFam;
    std::string prefix;

    bool snpMajorMode;
};

#endif /* _PLINKINPUTFILE_H_ */
