#ifndef _PLINKOUTPUTFILE_H_
#define _PLINKOUTPUTFILE_H_

#include <string>
#include <vector>

class SimpleMatrix;
class PlinkInputFile;
class VCFHeader;
class VCFRecord;

/****************************/
/*    Binary PLINK format   */
/****************************/
/*
 * Documentation
 * http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
 * BED (binary PED):
 * BIM (extended MAP file): chromosome, SNP, cM, base-position, allele 1, allele
 2
 *    e.g.
 *    1       snp1    0       1       G       A
 * FAM (first 6 columns of PED file)
 *      Family ID
 Individual ID
 Paternal ID
 Maternal ID
 Sex (1=male; 2=female; other=unknown)
 Phenotype
 *    e.g.
 *    1 1 0 0 1 0
 *
 */
class PlinkOutputFile {
 public:
  PlinkOutputFile(const char* fnPrefix) { init(fnPrefix); }
  PlinkOutputFile(const std::string& fnPrefix) { init(fnPrefix.c_str()); }
  void init(const char* fnPrefix);
  ~PlinkOutputFile() { close(); }
  void close() {
    if (this->fpFam) {
      fclose(this->fpFam);
      this->fpFam = NULL;
    }
    if (this->fpBim) {
      fclose(this->fpBim);
      this->fpBim = NULL;
    }
    if (this->fpBed) {
      fclose(this->fpBed);
      this->fpBed = NULL;
    }
  }
  void writeHeader(const VCFHeader* h);
  // @pos is from 0 to 3
  void setGenotype(unsigned char* c, const int pos, const int geno) {
    (*c) |= (geno << (pos << 1));
  }

  int writeRecord(VCFRecord* r);

  int writeRecordWithFilter(VCFRecord* r, const double minGD,
                            const double minGQ);
  /**
   * @return 0: success
   */
  int writeBIM(const char* chr, const char* id, double mapDist, int pos,
               const char* ref, const char* alt);

  /**
   * @return 0: success
   */
  int writeBIM(const std::vector<std::string>& chr,
               const std::vector<std::string>& id,
               const std::vector<double>& mapDist, const std::vector<int>& pos,
               const std::vector<std::string>& ref,
               const std::vector<std::string>& alt);

  void writeFAM(const std::string& people);
  void writeFAM(const std::vector<std::string>& people);
  void writeFAM(const std::vector<std::string>& fid,
                const std::vector<std::string>& iid,
                std::vector<double>& pheno);
  void writeFAM(const PlinkInputFile& pin, int idx);

  // NOTE: m should be: marker x people
  void writeBED(SimpleMatrix* mat, int nPeople, int nMarker);

  /**
   * Extract a subset of SNP and/or samples for binary PLINK file with @param
   * prefix
   * @param snpIdx indices for markers (should be >= 0)
   * @param sampleIdx indices for samples (should be >=0)
   */
  int extract(const std::string& prefix, const std::vector<int>& sampleIdx,
              const std::vector<int>& snpIdx);
  int extract(PlinkInputFile& pin, const std::vector<int>& sampleIdx,
              const std::vector<int>& snpIdx);
  int extractFAM(PlinkInputFile& pin, const std::vector<int>& sampleIdx);
  int extractFAMWithPhenotype(PlinkInputFile& pin,
                              const std::vector<int>& sampleIdx,
                              const SimpleMatrix& pheno);
  int extractBIM(PlinkInputFile& pin, const std::vector<int>& sampleIdx);
  int extractBED(PlinkInputFile& pin, const std::vector<int>& sampleIdx,
                 const std::vector<int>& snpIdx);

 private:
  int isMultiAllelic(const char* r);

 private:
  // we reverse the two bits as defined in PLINK format,
  // so we can process 2-bit at a time.
  const static unsigned char HOM_REF = 0x0;  // 0b00 ;
  const static unsigned char HET = 0x2;      // 0b10 ;
  const static unsigned char HOM_ALT = 0x3;  // 0b11 ;
  const static unsigned char MISSING = 0x1;  // 0b01 ;

  FILE* fpBed;
  FILE* fpBim;
  FILE* fpFam;
};  // end PlinkOutputFile

#endif /* _PLINKOUTPUTFILE_H_ */
