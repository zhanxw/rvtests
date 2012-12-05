#ifndef _GENOMESEQUENCE_H_
#define _GENOMESEQUENCE_H_

#include "base/IO.h"

// .fai format
// contig, size, location, basesPerLine, bytesPerLine
// e.g.
// 1       249250621       52      60      61
// 2       243199373       253404903       60      61
class Faidx{
public:
  struct Info{
    int contigSize;
    int offset;
    int basePerLine;
    int bytePerLine;
  };
public:
  /// @return: # of contigs read, minus number means errors!
  /// @param: file name
  int loadFaidx(const char* fn){
    LineReader lr(fn);
    std::vector<std::string> fd;
    int lineNo = 0;
    while (lr.readLineBySep(&fd, "\t")){
      lineNo ++;
      if (fd.size() != 5) {
        fprintf(stderr, "Wrong format: %s...\n", fd[0].c_str());
        continue;
      };
      Faidx::Info info;
      if ( !str2int( fd[1], &info.contigSize) ) {
        fprintf(stderr, "Cannot convert to integer at line %d!\n", lineNo);
        continue;
      }
      if ( !str2int( fd[2], &info.offset) ) {
        fprintf(stderr, "Cannot convert to integer at line %d!\n", lineNo);
        continue;
      }
      if ( !str2int( fd[3], &info.basePerLine) ) {
        fprintf(stderr, "Cannot convert to integer at line %d!\n", lineNo);
        continue;
      }
      if ( !str2int( fd[4], &info.bytePerLine) ) {
        fprintf(stderr, "Cannot convert to integer at line %d!\n", lineNo);
        continue;
      }
      if (this->data.count(fd[0]) != 0) {
        fprintf(stderr, "Warning, duplicate contig name at line %d!", lineNo);
      }
      this->data[fd[0]] = info;
    };
    return this->data.size();
  };
  Faidx::Info* getInfo(const std::string& chr){
    if (data.count(chr) == 0) {
      return NULL;
    } else {
      return &(data[chr]);
    }
  };
  int size() const {
    return data.size();
  };
  long int getGenomeLength() const {
    long int l = 0;
    std::map<std::string, Faidx::Info>::const_iterator it;
    for (it = data.begin(); it != data.end(); it++){
      l += it->second.contigSize;
    }
    return l;
  }
private:
  std::map<std::string, Faidx::Info> data;
};

class Chromosome{
public:
  explicit Chromosome(FILE* faFile, Faidx::Info* info): fp(faFile), info(info) {
  };
  explicit Chromosome():fp(NULL), info(NULL){};
  Chromosome(const Chromosome& chrom) {
    this->fp = chrom.fp;
    this->info = chrom.info;
  };
public:
  /**
   * Use 0-based index
   */
  char operator[] (unsigned int offset) const {
    int lineNo = offset / info->basePerLine;
    int remainder = offset % info->basePerLine;
    unsigned int pos = info->offset + lineNo * info->bytePerLine + remainder;
    if ( fseek(fp, pos, SEEK_SET) ) {
      fprintf(stderr, "Cannot fseek() at position %d!\n", pos);
      return 'N';
    }
    char c;
    if (1 != fread(&c, sizeof(char), 1, this->fp)) {
      fprintf(stderr, "Cannot fread()!\n");
      return 'N';
    }
    return c;
  };
  int size() const{
    return info->contigSize;
  };
private:
  FILE* fp;
  Faidx::Info* info;
};

class GenomeSequence{
public:
GenomeSequence():fp(NULL){};
  virtual ~GenomeSequence(){
    if (this->fp) {
      fclose(this->fp);
    }
  };
private:
  //forbid copying
  GenomeSequence(const GenomeSequence& gs);
  GenomeSequence& operator= (const GenomeSequence& gs);
public:
  /**
   * @return true: if loads successful
   */
  bool open(const char* fileName){
    // load .fa
    this->fp = fopen(fileName, "r");
    if (!this->fp) {
      fprintf(stderr, "Cannot open file: %s!", fileName);
      return false;
    }

    // load .fai
    std::string faiName = fileName;
    faiName.append(".fai");
    if (this->faidx.loadFaidx(faiName.c_str()) < 0) {
      fprintf(stderr, "Cannot open fai file!");
      return false;
    }
    
    return true;
  };
  
  /**
   * @return total number of chromosome
   */
  int size() const {
    return this->faidx.size();
  };
  /**
   * @return total number of chromosome
   */
  long int getGenomeLength() const {
    return faidx.getGenomeLength();
  };

  Chromosome& getChromosome(const std::string& c){
    std::map<std::string, Chromosome>::iterator it = data.find(c);
    if (it == data.end()) {
      Chromosome chrom(fp, faidx.getInfo(c));
      data[c] = chrom;
      return data[c];
    } else {
      return (it->second);
    }
  };
  const Chromosome& operator[] (const std::string& c) {
    std::map<std::string, Chromosome>::const_iterator it = data.find(c);
    if (it == data.end()) {
      Chromosome chrom(fp, faidx.getInfo(c));
      data[c] = chrom;
      return data[c];
    } else {
      return (it->second);
    }
  }
  /**
   * @return reference bases from chr:beg-end (beg inclusive, end exclusive)
   */
  std::string getBase(const std::string& c, int beg, int end) {
    std::string ret;
    if (beg > end) return ret;
    ret.resize( end - beg);
    for (int i = beg; i < end; ++i ) {
      ret[i- beg] = (*this)[c][i];
    };
    return ret;
  };
  bool exists(const std::string& c){
    if (this->data.find(c) != this->data.end())
      return true;
    if (faidx.getInfo(c) != NULL)
      return true;
    return false;
  }
public:
  std::map<std::string, Chromosome> data;
  FILE* fp;
  Faidx faidx;
};

#endif /* _GENOMESEQUENCE_H_ */
