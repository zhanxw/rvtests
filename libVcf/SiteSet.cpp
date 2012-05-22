#include "SiteSet.h"
#include "IO.h"
#include "TypeConversion.h"

/**
 * Load column 1 as chromosome, column 2 as position.
 * @return number of sites loaded
 */
int SiteSet::loadSiteFile(const char* fileName){
  int n = 0;
  LineReader lr(fileName);
  std::vector<std::string> fd;
  int pos;
  while(lr.readLineBySep(&fd, " \t")){
    if (fd.size() < 2) continue;
    pos = atoi(fd[1]);
    this->loadSite(fd[0], pos);
    ++n;
  };
  return n;
};

/**
 * Load plink .bim file
 * positions are 1-based index
 */
int SiteSet::loadBimFile(const char* fileName){
  int n = 0;
  LineReader lr(fileName);
  std::vector<std::string> fd;
  int pos;
  while(lr.readLineBySep(&fd, " \t")){
    if (fd.size() < 4) continue;
    pos = atoi(fd[3]);
    this->loadSite(fd[0], pos);
    ++n;    
  };
  return n;
};

/**
 * Load column 1 as chromosome, column 3 as position.
 * NOTE: rod file use 0-based index
 * @return number of sites loaded
 * example:
 1 rs55998931 10491
 1 rs62636508 10518
 1 rs58108140 10582
 1 rs10218492 10827
 1 rs10218493 10903
*/
int SiteSet::loadRodFile(const char* fileName){
  int n = 0;
  LineReader lr(fileName);
  std::vector<std::string> fd;
  int pos;
  while(lr.readLineBySep(&fd, " \t")){
    if (fd.size() < 3) continue;
    pos = atoi(fd[2]) + 1;
    this->loadSite(fd[0], pos);
    ++n;
  };
  return n;
};

// NOTE:
// BED file is [beg, end)
// here we don't care it is 0-based or 1-based and it's up to user to decide
int SiteSet::loadBEDFile(const char* fileName){
  int n = 0;
  LineReader lr(fileName);
  std::vector<std::string> fd;
  int beg;
  int end;
  while(lr.readLineBySep(&fd, " \t")){
    if (fd.size() < 3) continue;
    beg = atoi(fd[1]);
    end = atoi(fd[2]);
    for (int i = beg ; i < end; i++){
      this->loadSite(fd[0], i);
      ++n;
    }
  };
  return n;
};

