#ifndef _BGENVARIANT_H_
#define _BGENVARIANT_H_

#include <stdint.h>  // uint32_t
#include <cassert>
#include <string>
#include <vector>

class FileWriter;

static std::vector<std::vector<int> >
    table;  // for quickly lookup genotypes in colex order

struct BGenVariant {
  BGenVariant(const uint32_t& n) : N(n){};
  const uint32_t& N;
  // variant id data
  std::string varid;
  std::string rsid;
  std::string chrom;
  uint32_t pos;
  uint16_t K;                        // max = 2^16 = 65536
  std::vector<std::string> alleles;  // K alleles

  // store parsed results
  uint8_t B;  // bits to represent probability
  std::vector<bool> missing;
  std::vector<uint8_t> ploidy;  // Z
  bool isPhased;
  // prob[ index[0] .. index[1]] is the probability for the first individual
  // prob [ index[N-1] .. index[N]] is the probability for the last individual
  std::vector<int> index;
  std::vector<float> prob;  // probability array

  void makeTable(int ploidy, int allele) const;
  void findGenotype(int idx, int ploidy, int allele,
                    std::vector<int>* geno) const;
  void printGT(int i, FileWriter* fp) const;
  void printGTMissingFromHaplotype(FileWriter* fp) const;
  void printGTMissingFromGenotype(FileWriter* fp) const;
  void printGTAllele1FromGenotype(int i, FileWriter* fp) const;
  void printGTAllele2FromGenotype(int i, FileWriter* fp) const;
  void printGTAlleleGeneralFromGenotype(int idx, FileWriter* fp) const;
  void printGTFromHaplotype(int ii, FileWriter* fp) const;

  /// Handle GP  //////////////////////////////////////////////////
  /// output genotype probability
  void printGP(int i, FileWriter* fp) const;
  void printGPMissing(int idx, FileWriter* fp) const;
  void printGPAllele1(int i, FileWriter* fp) const;
  void printGPAllele2(int i, FileWriter* fp) const;
  void printGPAlleleGeneral(int idx, FileWriter* fp) const;
  /// Handle HP  //////////////////////////////////////////////////
  // handle haplotype probability
  void printHP(int i, FileWriter* fp) const;
  void printHPMissing(int idx, FileWriter* fp) const;
  void printHPAlleleGeneral(int idx, FileWriter* fp) const;

  /// Handle dosage //////////////////////////////////////////////////
  void printDosage(int i, FileWriter* fp) const;
};

#endif /* _BGENVARIANT_H_ */
