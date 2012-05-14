#ifndef _RANGELIST_H_
#define _RANGELIST_H_

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cassert>

#include "TypeConversion.h"
#include "Utils.h"
#include "IO.h"

/**
   range is inclusive on both edges.
*/
struct PositionPair{
    int begin;
    int end;
PositionPair():
    begin(0), end(0)
        {};
PositionPair(int b, int e):
    begin(b), end(e)
        {};
    bool operator== (const PositionPair o){
        return (begin == o.begin && end == o.end);
    }
    bool operator!= (const PositionPair o){
        return (begin != o.begin || end != o.end);
    }
};

inline bool PositionPairCompare(const PositionPair& p1, const PositionPair& p2){
    if (p1.begin != p2.begin)
        return (p1.begin < p2.begin);
    return (p1.end < p2.end);
};

/**
 * Store arbitrary number of ranges: e.g. [chr1:1-100, chr2:2-300....]
 */
class RangeCollection{
public:
    RangeCollection():_size(0){};
    void addRange(const std::string& chr, int begin, int end) {
        // if chr not exists
        // add chr to the chrVector
        // std::string c = chopChr(chr);
        if (this->rangeMap.find(chr) == this->rangeMap.end()){
            chrVector.push_back(chr);
        }

        // add begin, end to that chr entry
        this->rangeMap[chr].push_back( PositionPair(begin, end) );

        this->_size ++;
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

        this->_size = 0;
        for(int i = 0; i < this->chrVector.size(); i++ ) 
            this->_size += this->rangeMap[this->chrVector[i]].size();
    }

    void obtainRange(const int index, std::string* chrom, unsigned int* beg, unsigned int* end) const {
        int t = index;
        int s;
        for(int i = 0; i < this->chrVector.size(); i++ ) {
            const std::vector<PositionPair>& v = this->rangeMap.find(chrVector[i])->second;
            s = v.size();
            if ( t < s) {
                (*chrom) = chrVector[i];
                (*beg) = v[t].begin;
                (*end) = v[t].end;
                return;
            } else {
                t -= s;
            }
        }
        assert(false);
    };

    void obtainRange(const int index, std::string* range) const {
        unsigned int beg, end;
        this->obtainRange(index, range, &beg, &end);
        range->push_back(':');
        (*range) += toString(beg);
        range->push_back('-');
        (*range) += toString(end);
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
    void clear() {
        this->chrVector.clear();
        this->rangeMap.clear();
        this->_size = 0;
    };
    /**
     * @return number of regions
     */
    size_t size() const { 
        return this->_size;
    };
    const std::vector<std::string>& getChromVector() const{
        return this->chrVector;
    }
    const std::map< std::string, std::vector<PositionPair> >& getRangeMap() const{
        return this->rangeMap;
    }
  private:
    struct CompareChromName {
        bool operator() (const std::string& chrom1, const std::string& chrom2) { 
            return ( chrom2int(chrom1)  < chrom2int(chrom2));
        }
    } compareChromName;

    void sortChrVector() {
        dumpStringVector(chrVector);
        std::sort(chrVector.begin(), chrVector.end(), compareChromName);
        dumpStringVector(chrVector);
    };
    void dump(const std::vector<PositionPair>& v){
        for (int i = 0; i < v.size(); i++){
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
        int l = v->size();
        if (l == 0)
            return;

        int beg =  (*v)[0].begin;
        t.push_back( (*v)[0] );
        for(int i = 1; i < l; i++) {
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
private:
    std::vector<std::string> chrVector;
    std::map< std::string, std::vector<PositionPair> > rangeMap;
    size_t _size;
}; // end class RangeCollection


/**
 * Hold a collection of ranges, and offer utility function (e.g. load file)
 */
class RangeList{
public:
RangeList(): isSorted(false) {};
    void sort() {
        if (!this->isSorted) {
            this->rangeCollection.sort();
            this->isSorted = true;
        }
    }
    /**
     * First sort all ranges(if necessary), then store th index-th range to @param range, in the format of "chr1:100-200"
     */
    void obtainRange(const unsigned int index, std::string* range) {
        this->rangeCollection.obtainRange(index, range);
    }
    void obtainRange(const unsigned int index, std::string* chr, unsigned int* beg, unsigned int* end) {
        this->rangeCollection.obtainRange(index, chr, beg, end);
    }
    /**
     * @return total number of regions (after merging).
     */
    unsigned int size() const {return this->rangeCollection.size(); };
    // read gene list file and add these ranges
    void filterGeneName(const char* geneName, const char* fileName);
    /// argRangeList is a string indicating the range
    void addRangeList(const char* argRangeList);
    void addRangeFile(const char* argRangeFile);
    void addRange(const std::string& chr, int begin, int end) {
        this->isSorted = false;
        this->rangeCollection.addRange(chr, begin, end);
    };
    bool isInRange(const std::string& chr, int pos) {
        return this->rangeCollection.isInRange(chr, pos);
    };
    void clear() {
        this->rangeCollection.clear();
        this->isSorted = false;
    };

    
    class iterator{
   public:
  iterator(const RangeCollection& rc, int chrIdx, int inChrIdx):
        rangeCollection(&rc),
                chromIndex(chrIdx),
                chromSize(0),
                inChromRangeIndex(inChrIdx),
                inChromRegionSize(0),
                chrom(NULL),
                positionPair(NULL){
            chromSize = rc.getChromVector().size();
            inChromRegionSize = getRegions(chromIndex)? getRegions(chromIndex)->size() : 0;
            chrom = &(rc.getChromVector()[chromIndex]);
            positionPair = getRegions(chromIndex);
        };
        iterator& operator++(){
            this->inChromRangeIndex ++;
            if (inChromRangeIndex == inChromRegionSize) {
                this->chromIndex ++;
                this->inChromRangeIndex = 0;
                chromSize = this->rangeCollection->getChromVector().size();
                inChromRegionSize = getRegions(chromIndex)? getRegions(chromIndex)->size() : 0;
                chrom = &(this->rangeCollection->getChromVector()[chromIndex]);
                positionPair = getRegions(chromIndex);
            }
        }
        bool operator==(const iterator& iter) const{
            return (this->rangeCollection == iter.rangeCollection && 
                    this->chromIndex == iter.chromIndex &&
                    this->inChromRangeIndex == iter.inChromRangeIndex);
        };
        bool operator!=(const iterator& iter) const{
            return (this->rangeCollection != iter.rangeCollection ||
                    this->chromIndex != iter.chromIndex ||
                    this->inChromRangeIndex != iter.inChromRangeIndex);
        };
        const std::string getChrom() const{
            return (*this->chrom);
        };
        const unsigned int getBegin() const {
            return (*this->positionPair)[inChromRangeIndex].begin;
        };
        const unsigned int getEnd() const {
            return (*this->positionPair)[inChromRangeIndex].end;
        }
      private:
        bool isValidChromIndex(int index) const{
            if (0 <= index  && index < this->rangeCollection->getChromVector().size()) {
                return true;
            }
            return false;
        };
        const std::vector<PositionPair>* getRegions(const std::string& chrom) const{
            std::map< std::string, std::vector<PositionPair> >::const_iterator iter;
            iter = this->rangeCollection->getRangeMap().find(chrom);
            if (iter != this->rangeCollection->getRangeMap().end()) {
                return &(iter->second);                
            }
            fprintf(stderr, "ERROR: Don't contain chromosome %s!\n", chrom.c_str());
            return NULL;
        };
        const std::vector<PositionPair>* getRegions(int index) const{
            if (!isValidChromIndex(index) && index != this->rangeCollection->getChromVector().size()) {
                fprintf(stderr, "Invalid chromosome index: %d!\n", index);
                return NULL;
            }
            if (index == this->rangeCollection->getChromVector().size()) {
                return NULL;
            }
            const std::string& c = this->rangeCollection->getChromVector()[index];
            return getRegions(c);
        };

      private:
        const RangeCollection* rangeCollection;
        int chromIndex;
        int chromSize;
        int inChromRegionSize;
        int inChromRangeIndex;
        const std::string* chrom;
        const std::vector<PositionPair>* positionPair;
    };

    iterator begin() {
        return iterator(rangeCollection, 0, 0);
    };
    iterator end() {
        return iterator(rangeCollection, this->rangeCollection.getChromVector().size(), 0);
    };

private:
    RangeCollection rangeCollection;
    bool isSorted;
};


#endif /* _RANGELIST_H_ */
