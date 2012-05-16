#ifndef _SITESET_H_
#define _SITESET_H_

#include <string>
#include <set>
#include <map>

class SiteSet{
public:
    int addSiteFile(const char* fileName);
    int addSiteFile(const std::string& fileName) {
        return this->addSiteFile(fileName.c_str());
    };
    int addBEDFile(const char* fileName);
    int addBEDFile(const std::string& fileName) {
        return this->addBEDFile(fileName.c_str());
    };

    void addSite(const char* chrom, int pos) {
        site[chrom].insert(pos);
    };
    void addSite(const std::string& chrom, int pos) {
        site[chrom].insert(pos);
    };
    bool isIncluded(const char* chrom, int pos){
        this->_chrom = chrom;
        this->_it = site.find(_chrom);
        if (_it == site.end()) {
            return false;
        }
        if (_it->second.find(pos) == _it->second.end()){
            return false;
        }
        return true;
    }
    void clear() {
        this->site.clear();
    }
private:
    std::string _chrom;    //mutatable obj, just for speed up querying.
    std::map<std::string, std::set<int> >::const_iterator _it;
    std::map<std::string, std::set<int> > site;
};

#endif /* _SITESET_H_ */
