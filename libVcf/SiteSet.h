#ifndef _SITESET_H_
#define _SITESET_H_

#include <string>
#include <set>
//#include <map>
#include <unordered_set>
#include <unordered_map>
#include <cmath>

class SiteSet{
public:
  SiteSet() {
    this->site.rehash(ceil(32/this->site.max_load_factor()));
  }
  // Load plain position file
  // column 1: chrom, column 2: pos
  int loadSiteFile(const char* fileName);
  int loadSiteFile(const std::string& fileName) {
    return this->loadSiteFile(fileName.c_str());
  };
  // Load GATK .rod files
  int loadRodFile(const char* fileName);
  int loadRodFile(const std::string& fileName) {
    return this->loadRodFile(fileName.c_str());
  }
  // Load Plink .bim files
  int loadBimFile(const char* fileName);
  int loadBimFile(const std::string& fileName) {
    return this->loadBimFile(fileName.c_str());
  }
  // Load BED file
  int loadBEDFile(const char* fileName);
  int loadBEDFile(const std::string& fileName) {
    return this->loadBEDFile(fileName.c_str());
  };

  void loadSite(const char* chrom, int pos) {
    if (site.find(chrom) == site.end()) {
      site[chrom].rehash( ceil(1000000/site[chrom].max_load_factor()) );
    }
    site[chrom].insert(pos);
  };
  void loadSite(const std::string& chrom, int pos) {
    site[chrom].insert(pos);
  };
  bool isIncluded(const char* chrom, int pos) const{
    std::unordered_map<std::string, std::unordered_set<int> >::const_iterator it;
    it = this->site.find(chrom);
    if (it == this->site.end()) {
      return false;
    }
    if (it->second.find(pos) == it->second.end()){
      return false;
    }
    return true;
  }
  void clear() {
    this->site.clear();
  }
  size_t getTotalSite() const {
    size_t s = 0;
    std::unordered_map<std::string, std::unordered_set<int> >::const_iterator it = this->site.begin();
    for (; it != this->site.end(); it++) {
      s += it->second.size();
    }
    return s;
  }
private:
  std::unordered_map<std::string, std::unordered_set<int> > site;
};

#endif /* _SITESET_H_ */
