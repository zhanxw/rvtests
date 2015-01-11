#ifndef _ORDEREDMAP_H_
#define _ORDEREDMAP_H_

#include <stdio.h>
#include <cassert>
#include <vector>
#include <map>
#include "Exception.h"

/**
 * a simple OrderedMap class
 * use operator[] to insert elements
 * use at() to access elements
 * NOTE: KEY will be stored twice.
 */
template <class KEY, class TYPE>
    class OrderedMap{
public:
  bool find(const KEY& key) const {
    if (this->keyTypeMap.find(key) == this->keyTypeMap.end()){
      return false;
    }
    return true;
  }
  TYPE& operator[] (const KEY& key) {
    if (!this->find(key)){
      this->keyVec.push_back(key);
    }
    return this->keyTypeMap[key];
  }
  const TYPE& operator[] (const KEY& key) const{
    if (!this->find(key)){
      throw "key not found in OrderedMap";
    }
    return this->keyTypeMap.find(key)->second;
  }
  void front(KEY* k, TYPE* v) {
    *k = this->keyVec.front();
    *v = this->keyTypeMap[(*k)];
  }
  bool at(unsigned int idx, KEY* k, TYPE* v) {
    if (idx >= this->size()) return false;
    *k = this->keyVec[idx];
    *v = this->keyTypeMap[(*k)];
    return true;
  }
  bool at(unsigned int idx, KEY* k, TYPE* v) const {
    if (idx >= this->size()) return false;
    *k = this->keyVec[idx];
    if (this->keyTypeMap.find(*k) == this->keyTypeMap.end()){
      v =NULL;
    } else {
      *v = this->keyTypeMap.find(*k)->second;
    }
    return true;
  }
  const KEY& keyAt(unsigned int idx) const {
    if (idx >= this->size()) {
      log_error("Index out of bound, now quitting...");
    }
    return this->keyVec[idx];
  }
  const TYPE& valueAt(unsigned int idx) const {
    if (idx >= this->size()) {
      fprintf(stderr, "Cannot find KEY in valueAt()\n");
      abort();
    }
    const KEY& k = this->keyVec[idx];
    if (this->keyTypeMap.find(k) == this->keyTypeMap.end()){
      fprintf(stderr, "Cannot find KEY in valueAt()\n");
      abort();
    } else {
      return this->keyTypeMap.find(k)->second;
    }
  }
  TYPE& valueAt(unsigned int idx) {
    if (idx >= this->size()) {
      fprintf(stderr, "Cannot find KEY in valueAt()\n");
      abort();
    }
    const KEY& k = this->keyVec[idx];
    if (this->keyTypeMap.find(k) == this->keyTypeMap.end()){
      fprintf(stderr, "Cannot find KEY in valueAt()\n");
      abort();
    }
    return this->keyTypeMap.find(k)->second;

  }
  /**
   * compare
   */
  void compareKey(const OrderedMap<KEY, TYPE>& other, int* overlap, int* thisUniqueKeys, int* otherUniqueKeys) const{
    assert(overlap && thisUniqueKeys && otherUniqueKeys);
    *overlap = *thisUniqueKeys = *otherUniqueKeys = 0;
    for (unsigned int i = 0; i != this->size(); i ++ ){
      KEY& k = this->keyAt(i);
      if (other.find(k))
        (*overlap)++;
      else
        (*thisUniqueKeys)++;
    }
    *otherUniqueKeys = other.size() - *overlap;
  }
  unsigned int size() const { return this->keyVec.size();} ;
  void clear() {
    this->keyVec.clear();
    this->keyTypeMap.clear();
  };
private:
  std::vector < KEY > keyVec;
  std::map < KEY, TYPE > keyTypeMap;
};

#endif /* _ORDEREDMAP_H_ */
