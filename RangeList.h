#ifndef _RANGELIST_H_
#define _RANGELIST_H_

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cassert>

#include "TypeConversion.h"
#include "Utils.h"

/**
   range is inclusive on both edges.
 */
struct PositionPair{
    unsigned int begin;
    unsigned int end;
PositionPair():
    begin(0), end(0) 
        {};
PositionPair(unsigned int b, unsigned int e):
    begin(b), end(e) 
        {};
    bool operator== (const PositionPair o){
        return (begin == o.begin && end == o.end);
    }
};
bool PositionPairCompare(const PositionPair& p1, const PositionPair& p2){
    if (p1.begin != p2.begin)
        return (p1.begin < p2.begin);
    return (p1.end < p2.end);
}


class RangeCollection{
  public:
    void addRange(const std::string& chr, unsigned int begin, unsigned int end) {
        // if chr not exists
        // add chr to the chrVector
        std::string c = chopChr(chr);
        if (this->rangeMap.find(c) == this->rangeMap.end()){
            chrVector.push_back(c);
        }
        
        // add begin, end to that chr entry
        this->rangeMap[c].push_back( PositionPair(begin, end) );
    }
    void sort() {
        /*
         * some test code
         */
        /*
        this->addRange("X", 2, 4 );
        this->addRange("X", 1, 3 );
        this->addRange("1", 1, 3 );
        this->addRange("1", 4, 6 );
        this->addRange("1", 5, 10 );
        */
        this->sortChrVector();
        this->consolidate();
        /*
        assert(chrVector[0] == "1");
        assert(rangeMap["1"][0] == PositionPair(1, 3));
        assert(rangeMap["1"][1] == PositionPair(4, 10));
        assert(rangeMap["X"][0] == PositionPair(1, 4));
        */
    }
    
    void obtainRange(const unsigned int index, std::string* range) {
        int t = index; 
        int s; 
        for(unsigned int i = 0; i < this->chrVector.size(); i++ ) {
            s = rangeMap[chrVector[i]].size();
            if ( t < s) {
                (*range) = chrVector[i];
                range->push_back(':');
                (*range) += toString((this->rangeMap[this->chrVector[i]][t]).begin);
                range->push_back('-');
                (*range) += toString((this->rangeMap[this->chrVector[i]][t]).end);
                return;
            } else {
                t -= s;
            }
        }
        assert(false);
    };   
    size_t size() const { return this->rangeMap.size();};
  private:
    void sortChrVector() {
        std::sort(chrVector.begin(), chrVector.end());
    };
    void dump(const std::vector<PositionPair>& v){
        for (unsigned int i = 0; i < v.size(); i++){
            printf("[%d, %d] ", v[i].begin, v[i].end);
        }
        printf("\n");
    }
    void consolidate(){
        std::map< std::string, std::vector<PositionPair> >::iterator iter;
        for (iter = this->rangeMap.begin(); iter != this->rangeMap.end(); iter++) {
            sortPositionRange( & ((*iter).second));
            // dump( ((*iter).second));
            consolidateRange(& ((*iter).second) );
            // dump( ((*iter).second));
        }
    };
    // sort range in ascending order
    void sortPositionRange(std::vector<PositionPair>* v) {
        std::sort(v->begin(), v->end(), PositionPairCompare);
    };
    // now the range in v should be ordered
    // we will merge overlapped ranges
    void consolidateRange(std::vector<PositionPair>* v) {
        std::vector<PositionPair> t;
        unsigned int l = v->size();
        if (l == 0) 
            return;

        unsigned int beg =  (*v)[0].begin;
        t.push_back( (*v)[0] );
        for(unsigned int i = 1; i < l; i++) {
            // if this range falls into last range, skip
            if ( (*v)[i].end <= (*v)[i - 1].end )
                continue;
            
            // if this range overlaps with last range, extend last range
            if ( (*v)[i].begin <= (*v)[i-1].end && 
                 (*v)[i].end > (*v)[i-1].end) {
                t[i-1].end = (*v)[i].end;
                continue;
            }

            // this range now is after the last range, so append it to the list
            t.push_back( (*v)[i] );
        }
        
        // copy t -> v
        v->clear();
        std::swap( t, *v);
    };
    bool isInRange(const std::string& chr, int pos){
        if (rangeMap.find(chr) == rangeMap.end()) return false;
        std::vector<PositionPair>& r = rangeMap[chr];
        PositionPair p(pos, pos);
        int low = lower_bound(r.begin(), r.end(), p, PositionPairCompare) - r.begin();
        int up = upper_bound(r.begin(), r.end(), p, PositionPairCompare) - r.begin();
        for (int i = low; i < up; i++){
            PositionPair& pp = r[i];
            if (pp.begin <= pos && pp.end <= pos) {
                return true;
            }
        };
        return false;
    };
  private:
    std::vector<std::string> chrVector;
    std::map< std::string, std::vector<PositionPair> > rangeMap;
};

class RangeList{
public:
  RangeList(): isSorted(false) {};
    void obtainRange(const unsigned int index, std::string* range) {
        if (!this->isSorted) {
            this->rangeCollection.sort();
            this->isSorted = true;
        }
        this->rangeCollection.obtainRange(index, range);
    }
    unsigned int size() const {return this->rangeCollection.size(); };
    // read gene list file and add these ranges
    void filterGeneName(const char* geneName, const char* fileName);
    void addRangeList(const char* argRangeList);
    void addRangeFile(const char* argRangeFile);
    void addRange(const std::string& chr, unsigned int begin, unsigned int end) {
        this->isSorted = false;
        this->rangeCollection.addRange(chr, begin, end);
    };
    bool isInRange(const std::string& chr, int pos);
private:
    RangeCollection rangeCollection;
    bool isSorted;
};

//
void RangeList::filterGeneName(const char* inclusionGeneFileName, const char* geneTableFileName){
    // require user input gene list file
    if (strlen(geneTableFileName) == 0 && strlen(inclusionGeneFileName) != 0) {
        fprintf(stderr, "Please provide gene list file (e.g. refFlat) until we are able to process gene\n");
        exit(1);
    }

    // if not specify any gene, return whole range.
    if (strlen(inclusionGeneFileName) == 0) {
        return;
    }

    // store which gene do we want if specified
    std::set< std::string > inclusionSet;
    LineReader lr(inclusionGeneFileName);
    std::string gene;
    while(lr.readLine(&gene)) {
        inclusionSet.insert(gene);
    }
    
    std::vector<std::string> fields;
    std::string chr;
    std::string geneNameTbl;

    LineReader geneTable(geneTableFileName);
    while (geneTable.readLineBySep(&fields, "\t ")) {
        geneNameTbl = fields[0];
        if (inclusionSet.find(geneNameTbl) != inclusionSet.end()){ // store gene range
            chr = chopChr(fields[2].c_str());
            this->rangeCollection.addRange(chr, 
                                     atoi(fields[4].c_str()),   // start
                                     atoi(fields[5].c_str()));   // end
        }
    }
    if (this->rangeCollection.size() == 0){
        fprintf(stdout, "We cannot find given gene in your geneListFile, so all sites will be outputed\n");
    }
}

// verify if s is of the format: chr:begin-end format
// @return true: valid format
bool verifyRangeFormat(const std::string& s, char* chr, unsigned int* begin, unsigned int* end) {
    char c[128]; 
    char b[128];
    char e[128];
        
    int n = sscanf(s.c_str(), "%s:%s-%s", c, b, e);
    if ( n != 3 ) return false;
    int temp;
    if (!str2int(b, &temp) || temp <0) return false;
    int temp2;
    if (!str2int(e, &temp2) || temp2 <0) return false;
    if (temp > temp2) return false;
    return true;
}

// input range such as: 1:100-200,3:200-300
void RangeList::addRangeList(const char* argRangeList) {
    if (!strlen(argRangeList)) return;

    std::string rangeList;
    std::vector<std::string> col;
    //col.AddTokens(arg, ',');
    stringTokenize(rangeList, ',', &col);
    for (int i = 0; i < col.size(); i++){
        char c[64];
        unsigned int b,e;
        if (!verifyRangeFormat(col[i], c, &b, &e)) {
            fprintf(stdout, "This range does not conform 1:100-200 format -- %s\n", col[i].c_str());
        } else {
            this->rangeCollection.addRange(c, b, e);
        }
    }
};

/**
 * read a range list file like following
 * chr beg start
 * or 
 * chr beg
 * we will assume beg == end in the second case
 */
void RangeList::addRangeFile(const char* argRangeFile){
    if (!strlen(argRangeFile)) return;

    LineReader lr(argRangeFile);
    std::vector<std::string> sa;
    while ( lr.readLineBySep(&sa, "\t ")) {
        if (sa.size() == 0) continue;
        if (sa.size() == 1){
            fprintf(stderr, "Wrong --rangeFile: %s\n", argRangeFile);
            return;
        } else if (sa.size() == 2)
            this->rangeCollection.addRange(sa[0].c_str(), (unsigned int) atoi(sa[1]), (unsigned int) atoi(sa[1]));
        else if (sa.size() == 3)
            this->rangeCollection.addRange(sa[0].c_str(), (unsigned int) atoi(sa[1]), (unsigned int) atoi(sa[2]));
        else {
            fprintf(stdout, "Will only use the first 3 column of --rangeFile %s\n", argRangeFile);
            this->rangeCollection.addRange(sa[0].c_str(), (unsigned int) atoi(sa[1]), (unsigned int) atoi(sa[2]));
        }
    }
};

/**
 * check if chr:pos is in the collection
 */
bool RangeList::isInRange(const std::string& chr, int pos){
    
};

#endif /* _RANGELIST_H_ */
