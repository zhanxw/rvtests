#ifndef _PARREGION_H_
#define _PARREGION_H_

#include <limits.h>
#include <assert.h>
#include <map>
#include <set>
#include <string>
#include <utility> // std::pair
#include <vector>
#include "Utils.h"
#include "TypeConversion.h"

/**
 * A tools class to check if some variant in ParRegion
 * Refer: http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/
 */
class ParRegion{
 public:
  ParRegion() {
    init("chrX,X,chr23,23", "hg19");
  }
  ParRegion(const std::string& xLabel,
            const std::string& parRegion) {
    init(xLabel, parRegion);
  }
  int init(const std::string& xLabel,
           const std::string& parRegion) {
    if (xLabel.empty()) {
      // use default PLINK setting
      this->label.insert("X");
      this->label.insert("23");
    } else {
      std::vector<std::string> fd;
      stringTokenize(xLabel, ",", &fd);
      this->label.insert(fd.begin(), fd.end());
    }

    std::string region = tolower(parRegion);
    if (region.empty()) {
      region = "hg19";
    }
    if (region == "hg19" || region == "b37" || region == "grch37") {
      this->region.push_back(std::make_pair(60001, 2699520));
      this->region.push_back(std::make_pair(154931044, 155270560));
    } else if (region == "hg18" || region == "b36" || region == "grch36") {
      this->region.push_back(std::make_pair(1, 2709520));
      this->region.push_back(std::make_pair(154584238, 154913754));
    } else {
      std::vector<std::string> fd;
      std::vector<std::string> loc;
      stringTokenize(parRegion, ",", &fd);
      for (size_t i = 0; i < fd.size(); ++i) {
        stringTokenize(fd[i], "-", &loc);
        if (loc.size() != 2) {
          continue;
        }
        if (loc[0].empty()){ // skip region like "-100"
          continue;
        }
        int beg = atoi(loc[0]);
        int end = atoi(loc[1]);
        if (loc[1].empty()){
          end = INT_MAX;
        }
        this->region.push_back(std::make_pair(beg,end));
      }
    }
    return 0;
  }
  bool isParRegion(const std::string& chromPos) const {
    size_t p = chromPos.find(':');
    if (p == std::string::npos) return false; // malformat chromPos
    std::string chrom = chromPos.substr(0, p);
    size_t p2 = chromPos.find(p+1, ':');
    if (p2 != std::string::npos) return false; // malformat chromPos
    int pos = atoi(chromPos.substr(p + 1));
    return isParRegion(chrom, pos);
  }
  /**
   * @return true if @param chrom @param pos is in PAR region
   */
  bool isParRegion(const std::string& chrom,
                   int pos) const {
    if (this->isParRegionChrom(chrom) &&
        this->inParRegion(pos))
      return true;
    return false;
  }
  bool isHemiRegion(const std::string& chrom,
                    int pos) const {
  if (this->isParRegionChrom(chrom) &&
      !this->inParRegion(pos))
      return true;
    return false;
  }
  bool isParRegionChrom(const std::string& chrom) const {
    return (label.count(chrom) > 0);
  }
  /**
   * @return PAR chrommosome names concatenated by ","
   */
  std::string getLabel() const{
    assert (!this->label.empty());
    std::set<std::string>::const_iterator it = label.begin();
    std::string s = *it;
    for (it = label.begin(); it != label.end(); ++it) {
      s += ",";
      s += *it;
    }
    return s;
  }
  /**
   * @return PAR regions concatenated by ","
   */
  std::string getRegion() const {
    std::string s;
    for (size_t i = 0; i < this->region.size(); ++i) {
      if (i) s += ",";
      s += toString(region[i].first);
      s += "-";
      s += toString(region[i].second);
    }
    return s;
  }
  bool inParRegion(int pos) const{
    for(size_t i = 0; i < region.size(); ++i) {
      if (pos >= this->region[i].first &&
          pos <= this->region[i].second)
        return true;
    }
    return false;
  }
 private:
  std::set<std::string> label;
  std::vector< std::pair<int, int> > region; // inclusive booundaries
};

#endif /* _PARREGION_H_ */
