#ifndef _PLINKINPUTFILE_H_
#define _PLINKINPUTFILE_H_

#include "PlinkConst.h"

class PlinkInputFile{
public:
    PlinkInputFile(const char* fnPrefix) {
        std::string prefix = fnPrefix;
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
        fread(&c, sizeof(char), 1, this->fpBed);
        if (c != magic1) {
            fprintf(stderr, "Magic number of binary PLINK file does not match!\n");
            abort();
        }
        int mageic2 = 0x1b; // 0b00011011;
        fread(&c, sizeof(char), 1, this->fpBed);
        if (c != magic2) {
            fprintf(stderr, "Magic number of binary PLINK file does not match!\n");
            abort();
        }

        // snp major mode
        const ROW_MODE = 0x01; //0b00000001;
        const COL_MODE = 0x00; 
        fread(&c, sizeof(char), 1, this->fpBed);
        if ( c == ROW_MODE) {
            this->rowMode = true;
        }else if ( c== COL_MODE) {
            this->rowMode = false;
        } else{
            fprintf(stderr, "Unrecognized major mode in binary PLINK file.\n");
            abort();
        }
    };
    ~PlinkOutputFile() {
        fclose(this->fpBed);
        fclose(this->fpBim);
        fclose(this->fpFam);
    };
    
    int readIntoMatrix(Matrix* m) {
        assert(m);

    };
private:
    // we reverse the two bits as defined in PLINK format, 
    // so we can process 2-bit at a time.
    const static unsigned char HOM_REF = 0x0;     //0b00;
    const static unsigned char HET = 0x2;         //0b10;
    const static unsigned char HOM_ALT = 0x3;     //0b11;
    const static unsigned char MISSING = 0x1;     //0b01;

    FILE* fpBed;
    FILE* fpBim;
    FILE* fpFam;

    bool rowMode = false;
};

#endif /* _PLINKINPUTFILE_H_ */
