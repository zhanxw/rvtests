#ifndef _KGGINPUTFILE_H_
#define _KGGINPUTFILE_H_

#include <map>
#include <set>
#include <string>
#include <vector>

class BufferedReader;

class KGGInputFile {
 public:
  KGGInputFile(const std::string& fnPrefix);
  KGGInputFile(const std::string& fnPrefix, const std::string& fnSuffix);
  ~KGGInputFile();
  int init(const std::string& fnPrefix, const std::string& fnSuffix);

  // @return false if reached end
  bool readRecord();

  //////////////////////////////////////////////////
  // Sample inclusion/exclusion
  void includePeople(const std::string& s);
  void includePeople(const std::vector<std::string>& v);
  void includePeopleFromFile(const char* fn);
  void includeAllPeople();
  void excludePeople(const std::string& s);
  void excludePeople(const std::vector<std::string>& v);
  void excludePeopleFromFile(const char* fn);
  void excludeAllPeople();

  // No range related function
  int setSiteFile(const std::string& fn);

  int getGenotype(int indvIdx);
  void getAllele(int indvIdx, int* a1, int* a2);

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
  int getNumIndv() const { return this->indv.size(); }
  int getNumSample() const { return this->indv.size(); }
  int getNumEffectiveSample() const;
  int getNumMarker() const { return this->snp2Idx.size(); }
  const std::vector<std::string>& getIndv() const { return this->indv; }
  const std::vector<std::string>& getSampleName() const { return this->indv; }
  void getIncludedSampleName(std::vector<std::string>* p) const;
  int getEffectiveIndex(int idx) const;

  const std::vector<std::string>& getIID() const { return this->indv; }
  const std::vector<std::string>& getChrom() const { return this->chrom; }
  const std::vector<std::string>& getMarkerName() const { return this->snp; }
  const std::vector<double>& getMapDist() const { return this->mapDist; }
  const std::vector<int>& getPosition() const { return this->pos; }
  const std::vector<std::string>& getRef() const { return this->ref; }
  const std::vector<std::vector<std::string> >& getAlt() const {
    return this->alt;
  }
  const std::vector<int>& getSex() const { return this->sex; }
  const std::vector<double>& getPheno() const { return this->pheno; }

 public:
  std::vector<std::string> chrom;
  std::vector<std::string> snp;
  std::vector<double> mapDist;
  std::vector<int> pos;
  std::vector<std::string> ref;
  std::vector<std::vector<std::string> > alt;

  std::vector<std::string> indv;  /// people ids
  std::vector<int> sex;
  std::vector<double> pheno;

 private:
  void buildUnphasedTable(int numAllele);
  void buildPhasedTable(int numAllele);

  // sample inclusion/exclusion related
  void setPeopleMask(const std::string& s, bool b);
  void setPeopleMaskFromFile(const char* fn, bool b);
  void setRangeMode();
  // range list related
  void buildEffectiveIndex();
  void warnUnsupported(const char* tag);

 private:
  typedef struct TwoChar { unsigned char x[2]; } TwoChar;
  std::map<std::string, int> snp2Idx;
  std::map<std::string, int> pid2Idx;

  BufferedReader* fpKed;
  BufferedReader* fpKim;
  BufferedReader* fpKam;
  std::string prefix;

  int bits;                           // bits[0, 1, ... (bits-1)] is a block
  std::vector<unsigned char> buffer;  // raw data from ked file
  std::vector<unsigned char> data;    // genotype data
  int variantIdx;

  bool phased;
  std::map<int, std::map<char, TwoChar> > unphasedTable;
  std::map<int, std::map<char, TwoChar> > phasedTable;

  std::vector<bool> sampleMask;  // true means exclusion
  std::vector<int> effectiveIndex;
  // allow chromosomal sites
  std::set<std::string> allowedSite;
};

#endif /* _KGGINPUTFILE_H_ */
