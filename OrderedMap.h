#ifndef _ORDEREDMAP_H_
#define _ORDEREDMAP_H_

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
    unsigned int size() const { return this->keyVec.size();} ;
  private:
    std::vector < KEY > keyVec;
    std::map < KEY, TYPE > keyTypeMap;
};

#endif /* _ORDEREDMAP_H_ */
