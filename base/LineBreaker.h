#ifndef _LINEBREAKER_H_
#define _LINEBREAKER_H_

#include <string>
#include <vector>

/**
 * A simple line breaker
 */
class LineBreaker{
public:
LineBreaker(int width):
  separator(" \t"), width(width) {
  };
  void setSeparator(const std::string& s) {
    this->separator = s;
  };
  void setContent(const std::string& s) {
    this->formatted.clear();
    // use minimal lenth algorithm
    // see: http://en.wikipedia.org/wiki/Word_wrap
    int spaceLeft = width;
    std::string word;
    std::string space = "";
    size_t pos = 0;
    std::string line;
    while (getWord(s, pos, &word, &space)) {
        int taken = word.size() + space.size();
      if (taken > spaceLeft) { // add line break
        formatted.push_back(line);
        line.clear();
        spaceLeft = width - (word.size());
        line += word;
      } else { // continue fill in line
        line += space;
        line += word;
        spaceLeft -= space.size() + word.size();
      }
      pos += word.size() + space.size();
      // getWord(s, pos, &word, &space);
    }
    if (!line.empty()) {
      formatted.push_back(line);
    }
  };
  bool getWord(const std::string& s,
               size_t pos,
               std::string* word,
               std::string* space) {
    word->clear();
    space->clear();
    if (pos == s.size() || pos == std::string::npos)
      return false;
    // check space
    size_t beg = pos;
    while (separator.find(s[beg]) != std::string::npos && beg < s.size()) {
      space->push_back(s[beg]);
      beg++;
    }
    if (beg == s.size())
      return false;
    while (separator.find(s[beg]) == std::string::npos && beg < s.size()) {
      word->push_back(s[beg]);
      beg++;
    }
    return true;
  };
  void clear() {
    this->formatted.clear();
  };
  size_t getWidth() const {
    return width;
  };
  size_t getHeight() const {
    return formatted.size();
  };
  size_t size() const {
    return formatted.size();
  };
  const std::string& operator[] (const int idx) const {
    return formatted[idx];    
  }
  const std::string& getLine(int idx)  const{
    return formatted[idx];
  };
private:
  std::string separator;
  int width;
  std::vector<std::string> formatted;
}; //end LineBreaker

void printTwoColumn(FILE* fp,  LineBreaker& left, LineBreaker& right, const char* sep) {
  size_t ls = left.size();
  size_t rs = right.size();
  size_t lmin = ls < rs ? ls : rs;
  for (size_t i = 0; i < lmin; ++i) {
    fprintf(fp, "%*s%s%s\n", (int)left.getWidth(), left[i].c_str(),
            sep,
            right[i].c_str());
  }
  if (ls > rs) {
    for (size_t i = lmin; i < ls; ++i) {
      fprintf(fp, "%s%s\n", left[i].c_str(), sep);
    }
  } else if (ls < rs){
    for (size_t i = lmin; i < rs; ++i) {
      fprintf(fp, "%*s%s%s\n", (int)left.getWidth(), "", sep, right[i].c_str());
    }
  }
};
#endif /* _LINEBREAKER_H_ */
