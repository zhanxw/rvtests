#ifndef _PLINKINPUTFILE_H_
#define _PLINKINPUTFILE_H_

#include <map>
#include <string>

#include "base/Exception.h"
#include "base/IO.h"

class SimpleMatrix;

class PlinkInputFile {
 public:
  PlinkInputFile(const std::string& fnPrefix) {
    this->prefix = fnPrefix;
    this->fpBed = fopen((prefix + ".bed").c_str(), "rb");
    this->fpBim = fopen((prefix + ".bim").c_str(), "rt");
    this->fpFam = fopen((prefix + ".fam").c_str(), "rt");
    if (!this->fpBed || !this->fpBim || !this->fpFam) {
      REPORT("Cannot open binary PLINK file!");
      abort();
    }
    // write Bed header
    char c;
    // magic number
    char magic1 = 0x6c;  // 0b01101100;
    int ret;
    UNUSED(ret);
    ret = fread(&c, sizeof(char), 1, this->fpBed);
    if (ret != 1) {
      fprintf(stderr, "Encounter error when reading plink BED files.\n");
    }
    assert(ret == 1);
    if (c != magic1) {
      fprintf(stderr, "Magic number of binary PLINK file does not match!\n");
      abort();
    }
    int magic2 = 0x1b;  // 0b00011011;
    ret = fread(&c, sizeof(char), 1, this->fpBed);
    if (ret != 1) {
      fprintf(stderr, "Encounter error when reading plink magic number.\n");
    }
    assert(ret == 1);
    if (c != magic2) {
      fprintf(stderr, "Magic number of binary PLINK file does not match!\n");
      abort();
    }

    // snp major mode
    const int SNP_MAJOR_MODE = 0x01;  // 0b00000001;
    const int INDV_MAJOR_MODE = 0x00;
    ret = fread(&c, sizeof(char), 1, this->fpBed);
    if (ret != 1) {
      fprintf(stderr, "Encounter error when determining plink files modes.\n");
    }
    assert(ret == 1);
    if (c == SNP_MAJOR_MODE) {
      this->snpMajorMode = true;
    } else if (c == INDV_MAJOR_MODE) {
      this->snpMajorMode = false;
    } else {
      fprintf(stderr, "Unrecognized major mode in binary PLINK file.\n");
      exit(1);
    }

    // read bim
    LineReader* lr = new LineReader((this->prefix + ".bim").c_str());
    std::string chrPos;
    std::vector<std::string> fd;
    while (lr->readLineBySep(&fd, " \t")) {
      if (fd.size() != 6) {
        fprintf(stderr, "Wrong format in bim file.\n");
        exit(1);
      }

      // when rsid == ".", use "chr:pos" as dict key
      if (fd[1] == ".") {
        chrPos = fd[0];
        chrPos += ":";
        chrPos += fd[3];
      } else {
        chrPos = fd[1];
      }
      if (snp2Idx.find(chrPos) == snp2Idx.end()) {
        chrom.push_back(fd[0]);
        snp.push_back(fd[1]);
        snp2Idx[chrPos] = 0;
        const int val = snp2Idx.size() - 1;
        snp2Idx[chrPos] = val;
        mapDist.push_back(atof(fd[2].c_str()));
        pos.push_back(atoi(fd[3].c_str()));
        ref.push_back(fd[4]);
        alt.push_back(fd[5]);
      } else {
        fprintf(stderr,
                "Error found: duplicated marker name or chromosomal position [ "
                "%s ]!\n",
                fd[1].c_str());
        exit(1);
      }
    }
    delete lr;

    // read fam
    lr = new LineReader((this->prefix + ".fam").c_str());
    while (lr->readLineBySep(&fd, " \t")) {
      if (fd.size() != 6) {
        fprintf(stderr, "Wrong format in fam file.\n");
        exit(1);
      }

      // will skip loading fam, fatherid, motherid
      if (pid2Idx.find(fd[1]) == pid2Idx.end()) {
        pid2Idx[fd[1]] = -1;
        const int val = pid2Idx.size() - 1;
        pid2Idx[fd[1]] = val;
        indv.push_back(fd[1]);
        sex.push_back(atoi(fd[4].c_str()));
        pheno.push_back(atof(fd[5].c_str()));
      } else {
        fprintf(stderr, "duplicated person id [ %s ], ignore!\n",
                fd[1].c_str());
        exit(1);
      }
    }
    delete lr;

    fprintf(stderr,
            "Finished loading %s.{fam,bim,bed}, %zu markers, %zu samples\n",
            fnPrefix.c_str(), snp2Idx.size(), indv.size());
  }
  ~PlinkInputFile() {
    fclose(this->fpBed);
    fclose(this->fpBim);
    fclose(this->fpFam);
  }

  // @param m: people by marker matrix
  int readIntoMatrix(SimpleMatrix* mat) const;
  int readIntoMatrix(SimpleMatrix* mat, std::vector<std::string>* peopleNames,
                     std::vector<std::string>* markerNames) const;

  // summary statistics calculations
  int calculateMAF(std::vector<double>* maf);
  int calculateMissing(std::vector<double>* imiss, std::vector<double>* lmiss);

  // read BED file
  int readBED(unsigned char* buf, size_t n);

  // utility functions
  // get PLINK 2bit genotype for the @param sample'th sample and @param
  // marker'th marker
  unsigned char get2BitGenotype(int sample, int marker);
  // @param m is the maker name.
  // can be used to check whether a marker exists
  int getMarkerIdx(const std::string& m) {
    if (this->snp2Idx.find(m) == this->snp2Idx.end()) {
      return -1;
    } else {
      return (this->snp2Idx[m]);
    }
  }
  // @param p is the sample name
  // can be used to check whether a sample exists
  int getSampleIdx(const std::string& p) {
    if (this->pid2Idx.find(p) == this->pid2Idx.end()) {
      return -1;
    } else {
      return (this->pid2Idx[p]);
    }
  }
  // get the @param i-th 2-bit genotype
  // NOTE: when g = |ddcc|bbaa| (8bits)
  // extract2Bit(g, 0) = |0000|00aa| (8bits)
  static unsigned char extract2Bit(unsigned char g, int i);
  int getNumIndv() const { return this->indv.size(); }
  int getNumSample() const { return this->indv.size(); }
  int getNumMarker() const { return this->snp2Idx.size(); }
  const std::vector<std::string>& getIndv() const { return this->indv; }
  const std::vector<std::string>& getSampleName() const { return this->indv; }
  const std::vector<std::string>& getIID() const { return this->indv; }
  const std::vector<std::string>& getChrom() const { return this->chrom; }
  const std::vector<std::string>& getMarkerName() const { return this->snp; }
  const std::vector<double>& getMapDist() const { return this->mapDist; }
  const std::vector<int>& getPosition() const { return this->pos; }
  const std::vector<std::string>& getRef() const { return this->ref; }
  const std::vector<std::string>& getAlt() const { return this->alt; }
  const std::vector<int>& getSex() const { return this->sex; }
  const std::vector<double>& getPheno() const { return this->pheno; }

 public:
  std::vector<std::string> chrom;
  std::vector<std::string> snp;
  std::vector<double> mapDist;
  std::vector<int> pos;
  std::vector<std::string> ref;
  std::vector<std::string> alt;

  std::vector<std::string> indv;  /// people ids
  std::vector<int> sex;
  std::vector<double> pheno;

 public:
  // we reverse the two bits as defined in PLINK format,
  // so we can process 2-bit at a time.
  const static unsigned char HOM_REF = 0x0;  // 0b00;
  const static unsigned char HET = 0x2;      // 0b10;
  const static unsigned char HOM_ALT = 0x3;  // 0b11;
  const static unsigned char MISSING = 0x1;  // 0b01;

 private:
  std::map<std::string, int> snp2Idx;
  std::map<std::string, int> pid2Idx;

  FILE* fpBed;
  FILE* fpBim;
  FILE* fpFam;
  std::string prefix;

  bool snpMajorMode;
};

#endif /* _PLINKINPUTFILE_H_ */
