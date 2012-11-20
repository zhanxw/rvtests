#include "RangeList.h"
#include <climits>

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

/**
 * verify if s is of the format: chr:begin-end format
 * @return 0: valid format
 */
int parseRangeFormat(const std::string& s, std::string* chr, unsigned int* begin, unsigned int* end) {
  //fprintf(stderr, "parseRangeFormat: %s\n", s.c_str());

  unsigned int i = 0;
  chr->clear();
  while (i < s.size()){
    if (s[i]!=':'){
      chr->push_back(s[i]);
    } else{
      break;
    }
    i++;
  }
  i ++; //skip ':'

  std::string t;
  while (i < s.size()){
    if (s[i] !='-') {
      t.push_back(s[i]);
    } else{
      break;
    }
    i++;
  }
  int b = 0;
  if (!str2int(t.c_str(), &b) || b < 0) return -1;
  *begin = b;

  if (s[i] == '\0'){ // 1:100 meaning from 1:100- 1:(max_pos)
    *end = 1 << 29;  /// this is copied from tabix index.c ti_parse_region
    // fprintf(stderr, "UINT_MAX = %u", UINT_MAX);
    return 0;
  }

  i ++ ; // skip '-'
  int e = 1<<29;     // that's the constant used in tabix
  if (s[i] == '\0'){ //format like: 1:100-
    *end = 1<<29;  // that's the constant used in tabix
  } else{
    if (!str2int(s.c_str() + i, &e) || e < 0 || b > e) return false;
  }
  *end = e;

  // fprintf(stderr, "parse result: %s %d %d\n", chr->c_str(), *begin, *end);
  return 0;
}

/**
 * input range such as:
 * 1:100-200,3:200-300
 * X:150
 * MT
 */
void RangeList::addRangeList(const char* argRangeList) {
  if (!strlen(argRangeList)) return;

  std::string rangeList = argRangeList;
  std::vector<std::string> col;
  //col.AddTokens(arg, ',');
  stringTokenize(rangeList, ',', &col);
  for (unsigned int i = 0; i < col.size(); i++){
    std::string c;
    unsigned int b,e;
    if (!parseRangeFormat(col[i], &c, &b, &e)) {
      this->rangeCollection.addRange(c, b, e);
    } else {
      fprintf(stdout, "This range does not conform 1:100-200 format -- skip %s\n", col[i].c_str());
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
  // fprintf(stdout, "Load range file %s.\n", argRangeFile);

  LineReader lr(argRangeFile);
  std::vector<std::string> sa;
  while ( lr.readLineBySep(&sa, "\t ")) {
    if (sa.size() == 0) continue;
    if (sa.size() == 1){
      // fprintf(stderr, "Wrong format for --rangeFile: %s, shoudl be: chr beg end \n", argRangeFile);
      this->addRangeList(sa[0].c_str());
    } else if (sa.size() == 2)
      this->rangeCollection.addRange(sa[0].c_str(), (unsigned int) atoi(sa[1]), (unsigned int) atoi(sa[1]));
    else if (sa.size() == 3)
      this->rangeCollection.addRange(sa[0].c_str(), (unsigned int) atoi(sa[1]), (unsigned int) atoi(sa[2]));
    else {
      // we will silently use the first 3 columns
      // fprintf(stdout, "Will only use the first 3 column of --rangeFile %s\n", argRangeFile);
      this->rangeCollection.addRange(sa[0].c_str(), (unsigned int) atoi(sa[1]), (unsigned int) atoi(sa[2]));
    }
  }
};

