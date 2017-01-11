#ifndef _INDEXMAP_H_
#define _INDEXMAP_H_

/**
 * IndexMap use integer index as keys (e.g. 0, 1, 2, ...) and customized type as
 * values
 */
template <class TYPE>
class IndexMap {
 public:
  TYPE& operator[](const int& key) {
    if (key < 0) {
      fprintf(stderr, "%s:%d Wrong key [ %d ]\n", __FILE__, __LINE__, key);
      exit(1);
    }
    if (key >= (int)value.size()) {
      value.resize(key + 1);
    }
    return value[key];
  }
  const TYPE& operator[](const int& key) const {
    if (key < 0 || key >= (int)value.size()) {
      fprintf(stderr, "%s:%d Wrong key [ %d ]\n", __FILE__, __LINE__, key);
      exit(1);
    }
    return value[key];
  }
  size_t size() const { return this->value.size(); }
  void clear() { this->value.clear(); }

 private:
  std::vector<TYPE> value;
};

#endif /* _INDEXMAP_H_ */
