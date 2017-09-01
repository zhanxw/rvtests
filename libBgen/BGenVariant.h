#ifndef _BGENVARIANT_H_
#define _BGENVARIANT_H_

#include <stdint.h>  // uint32_t
#include <cassert>
#include <string>
#include <vector>

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
  std::vector<bool> missing;
  std::vector<uint8_t> ploidy;  // Z
  bool isPhased;
  // prob[ index[0] .. index[1]] is the probability for the first individual
  // prob [ index[N-1] .. index[N]] is the probability for the last individual
  std::vector<int> index;
  std::vector<float> prob;  // probability array

  void makeTable(int ploidy, int allele) const {
    assert(ploidy >= 0 && allele >= 0);
    if ((size_t)ploidy > table.size()) {
      table.resize(ploidy);
    }
    for (int i = 0; i < ploidy; ++i) {
      const int s = table[i].size();
      if (allele > s) {
        table[i].resize(allele);
      }
      for (int j = s; j < allele; ++j) {
        if (i == 0) {
          table[i][j] = 1;
          continue;
        }

        if (j == 0) {
          table[i][j] = 1;
        } else if (j == 1) {
          table[i][j] = i + 1;
        } else {
          table[i][j] = table[i][j - 1] * (i + j) / j;
        }
      }
    }
  }
  void findGenotype(int idx, int ploidy, int allele,
                    std::vector<int>* geno) const {
    std::vector<int>& g = *geno;
    g.resize(ploidy);

    makeTable(ploidy, allele);
    --ploidy;
    ++idx;
    while (ploidy >= 0) {
      int i = 0;
      int cumsum = 0;
      for (; i < allele; ++i) {
        if (cumsum + table[ploidy][i] >= idx) {
          idx -= cumsum;
          break;
        } else {
          cumsum += table[ploidy][i];
        }
      }

      g[ploidy] = i;
      ploidy--;
    }
  }
  void printGT(int i, FILE* fp) const {
    if (isPhased) {
      if (missing[i]) {
        printGTMissingFromHaplotype(fp);
      } else {
        printGTFromHaplotype(i, fp);
      }
    } else {
      if (missing[i]) {
        printGTMissingFromGenotype(fp);
      } else {
        switch (K) {
          case 2:
            printGTAllele2FromGenotype(i, fp);
            break;
          case 1:
            printGTAllele1FromGenotype(i, fp);
            break;
          default:
            printGTAlleleGeneralFromGenotype(i, fp);
            break;
        }
      }
    }
  }
  void printGTMissingFromHaplotype(FILE* fp) const {
    fputs(".", fp);
    for (int i = 1; i < ploidy[i]; ++i) {
      fputs("|.", fp);
    }
  }
  void printGTMissingFromGenotype(FILE* fp) const {
    fputs(".", fp);
    for (int i = 1; i < ploidy[i]; ++i) {
      fputs("/.", fp);
    }
  }
  void printGTAllele1FromGenotype(int i, FILE* fp) const {
    fputs("0", fp);
    for (int i = 1; i < ploidy[i]; ++i) {
      fputs("/0", fp);
    }
  }
  void printGTAllele2FromGenotype(int i, FILE* fp) const {
    if (ploidy[i] == 2) {  // prob[index[i]] stores p(00), p(01), p(11)
      const float prob0 = prob[index[i]];
      const float prob1 = prob[index[i] + 1];
      const float prob2 = prob[index[i] + 2];

      if (prob0 > prob1 && prob0 > prob2) {
        fputs("0/0", fp);
      } else if (prob1 > prob0 && prob1 > prob2) {
        fputs("0/1", fp);
      } else {
        fputs("1/1", fp);
      }
    } else if (ploidy[i] == 1) {  // prob[index[i]] stores p(0), p(1)
      const float prob0 = prob[index[i]];
      const float prob1 = prob[index[i] + 1];

      if (prob0 > prob1) {
        fputs("0", fp);
      } else {
        fputs("1", fp);
      }
    } else {  // prob[index[i]] stores p(00..0), p(00...1), ... p(22..2)
      printGTAlleleGeneralFromGenotype(i, fp);
    }
  }
  void printGTAlleleGeneralFromGenotype(int idx, FILE* fp) const {
    // let Z = ploidity, K = alleles
    // total choose(N+K-1, K-1) genotypes
    int maxIdx = index[idx];
    float maxVal = prob[maxIdx];
    for (int i = maxIdx + 1; i < index[idx + 1]; ++i) {
      if (prob[i] > maxVal) {
        maxIdx = i;
        maxVal = prob[i];
      }
    }
    std::vector<int> geno;
    findGenotype(maxIdx - index[idx], ploidy[idx], K, &geno);
    for (size_t i = 0; i < geno.size(); ++i) {
      if (i) {
        fputs("/", fp);
      }
      fprintf(fp, "%d", geno[i]);
    }
  }
  void printGTFromHaplotype(int ii, FILE* fp) const {
    const int Z = ploidy[ii];
    //  K allels
    int idx = index[ii];
    int maxIdx;
    float maxValue;
    for (int i = 0; i < Z; ++i) {
      maxIdx = 0;
      maxValue = prob[idx];
      idx++;
      for (int j = 1; j < K; ++j) {
        if (prob[idx] > maxValue) {
          maxValue = prob[idx];
          maxIdx = j;
        }
        idx++;
      }
      if (i) {
        fputs("|", fp);
      }
      fprintf(fp, "%d", maxIdx);
    }
    assert(idx == index[ii + 1]);
  }

  /// Handle GP  //////////////////////////////////////////////////
  /// output genotype probability
  void printGP(int i, FILE* fp) const {
    if (missing[i]) {
      printGPMissing(i, fp);
      return;
    }

    switch (K) {
      case 2:
        printGPAllele2(i, fp);
        break;
      case 1:
        printGPAllele1(i, fp);
        break;
      default:
        printGPAlleleGeneral(i, fp);
        break;
    }
  }
  void printGPMissing(int idx, FILE* fp) const {
    for (int i = index[idx]; i < index[idx + 1]; ++i) {
      if (i != index[idx]) {
        fputs(",", fp);
      }
      fputs(".", fp);
    }
  }
  void printGPAllele1(int i, FILE* fp) const { fputs("1", fp); }
  void printGPAllele2(int i, FILE* fp) const {
    if (ploidy[i] == 2) {  // prob of 00, 01, 11
      const float prob0 = prob[index[i]];
      const float prob1 = prob[index[i] + 1];
      const float prob2 = prob[index[i] + 2];

      fprintf(fp, "%g,%g,%g", prob0, prob1, prob2);
    } else if (ploidy[i] == 1) {  // prob of 0,1
      const float prob0 = prob[index[i]];
      const float prob1 = prob[index[i] + 1];

      fprintf(fp, "%g,%g", prob0, prob1);
    } else {
      // let Z = ploidity, K = alleles
      // total choose(N+K-1, K-1) genotypes
      printGPAlleleGeneral(i, fp);
    }
  }
  void printGPAlleleGeneral(int idx, FILE* fp) const {
    for (int i = index[idx]; i < index[idx + 1]; ++i) {
      if (i != index[idx]) {
        fputs(",", fp);
      }
      fprintf(fp, "%g", prob[i]);
    }
  }
  /// Handle HP  //////////////////////////////////////////////////
  // handle haplotype probability
  void printHP(int i, FILE* fp) const {
    if (missing[i]) {
      printHPMissing(i, fp);
      return;
    }
    printHPAlleleGeneral(i, fp);
  }
  void printHPMissing(int idx, FILE* fp) const {
    for (int i = index[idx]; i < index[idx + 1]; ++i) {
      if (i != index[idx]) {
        fputs(",", fp);
      }
      fputs(".", fp);
    }
  }
  void printHPAlleleGeneral(int idx, FILE* fp) const {
    for (int i = index[idx]; i < index[idx + 1]; ++i) {
      if (i != index[idx]) {
        fputs(",", fp);
      }
      fprintf(fp, "%g", prob[i]);
    }
  }

  /// Handle dosage //////////////////////////////////////////////////
  void printDosage(int i, FILE* fp) const {
    if (missing[i]) {
      fputs(".", fp);
      return;
    }
    if (ploidy[i] == 2 && K == 2) {
      // const float prob0 = prob[index[i]];
      const float prob1 = prob[index[i] + 1];
      const float prob2 = prob[index[i] + 2];

      fprintf(fp, "%g", prob1 + 2.0 * prob2);
    } else {
      fputs(".", fp);
    }
  }
};

#endif /* _BGENVARIANT_H_ */
