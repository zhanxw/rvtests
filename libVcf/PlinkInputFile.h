#ifndef _PLINKINPUTFILE_H_
#define _PLINKINPUTFILE_H_

#include <map>
#include <string>

#include "Exception.h"
#include "IO.h"

class SimpleMatrix;

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

    fprintf(stderr, "Finished loading %s. %zu chrom, %zu indv\n", fnPrefix, snp2Idx.size(), indv.size());
  };
  ~PlinkInputFile() {
    fclose(this->fpBed);
    fclose(this->fpBim);
    fclose(this->fpFam);
  };

  // @param m: people by marker matrix
  int readIntoMatrix(SimpleMatrix* mat);
  int readIntoMatrix(SimpleMatrix* mat, std::vector<std::string>* peopleNames, std::vector<std::string>* markerNames);

  int getMarkerIdx(const std::string& m){
    if (this->snp2Idx.find(m) == this->snp2Idx.end()){
      return -1;
    } else{
      return (this->snp2Idx[m]);
    }
  };

  int getNumIndv() const {
    return this->indv.size();
  };
  int getNumMarker() const {
    return this->snp2Idx.size();
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
