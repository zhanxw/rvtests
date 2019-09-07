#include <stdlib.h>
#include <string>
#include <vector>

#include "base/Utils.h"
#include "third/samtools/bgzf.h"

// 1-st level keys: contig, FILTER, INFO, FORMAT
// 2-st level keys:
//   contig: ID, length
//   FILTER: ID, description, (will manually add Number=1, Type=String)
//   INFO: ID, Number, Type, Description
//   FORMAT: ID, Number, Type, Description
// we use a vector to store contigs (header_contig_id)
// and four vectors to store FILTER, INFO, FORMATs (header_id, header_number,
// header_type, header_description)
std::vector<std::string> header_contig_id;
std::vector<std::string> header_id;
std::vector<std::string> header_number;
std::vector<std::string> header_type;
std::vector<std::string> header_description;

class BCFHeaderParser {
 public:
  std::string parseValue(const std::string& s, const std::string& key) {
    size_t begin = s.find("<");
    size_t end = s.rfind(">");
    if (begin == std::string::npos || end == std::string::npos) {
      fprintf(stderr, "Wrong intput string during parsing!\n");
    }
    std::string ss = s.substr(begin + 1, end - begin - 1);
    // printf("ss = %s\n", ss.c_str());
    begin = ss.find(key);
    ss = ss.substr(begin, ss.size() - begin);
    if (ss.substr(0, key.size()) != key || ss[key.size()] != '=') {
      fprintf(stderr, "Cannot find the key\n");
    }
    ss = ss.substr(key.size() + 1, ss.size() - key.size() - 1);
    // printf("ss = %s\n", ss.c_str());
    if (ss[0] == '"') {
      end = ss.find("\"", 1);
      ss = ss.substr(1, end - 1);
    } else {
      end = ss.find_first_of(",>");
      ss = ss.substr(0, end);
    }

    // printf("ss = %s\n", ss.c_str());
    return ss;
  }
  int parse(const std::string& s) {
    key = id = number = type = desc = "";
    if (s.find("##contig=<") == 0) {
      key = "contig";
      id = parseValue(s, "ID");
    } else if (s.find("##FILTER=<") == 0) {
      key = "filter";
      id = parseValue(s, "ID");
      number = "1";
      type = "String";
      desc = parseValue(s, "Description");
      idx = atoi(parseValue(s, "IDX"));
    } else if (s.find("##INFO=<") == 0) {
      key = "info";
      id = parseValue(s, "ID");
      number = parseValue(s, "Number");
      type = parseValue(s, "Type");
      desc = parseValue(s, "Description");
      idx = atoi(parseValue(s, "IDX"));
    } else if (s.find("##FORMAT=<") == 0) {
      key = "format";
      id = parseValue(s, "ID");
      number = parseValue(s, "Number");
      type = parseValue(s, "Type");
      desc = parseValue(s, "Description");
      idx = atoi(parseValue(s, "IDX"));
    }
    return 0;
  }
  std::string key;
  std::string id;
  std::string number;
  std::string type;
  std::string desc;
  int idx;
};

int parseHeader(const std::string& header,
                std::vector<std::string>* header_contig_id,
                std::vector<std::string>* header_id,
                std::vector<std::string>* header_number,
                std::vector<std::string>* header_type,
                std::vector<std::string>* header_desc) {
  std::vector<std::string> lines;
  stringTokenize(header, "\n", &lines);
  BCFHeaderParser parser;
  for (size_t i = 0; i != lines.size(); ++i) {
    if (parser.parse(lines[i]) < 0) {
      fprintf(stderr, "Parser encountered error!\n");
      return -1;
    }
    if (parser.key == "contig") {
      header_contig_id->push_back(parser.id);
    } else if (parser.key == "filter" || parser.key == "info" ||
               parser.key == "format") {
      // need to check the hidden IDX= value
      // e.g. "ID=PASS,Description="All filters passed",IDX=0"
      // in BCF implementation, the default and hidden IDX=0 specified the
      // dictionary index
      // this is important as FILTEr and INFO can share the same ID, for
      // example:
      // ##FILTER=<ID=SVM,Description="Variant failed SVM filter",IDX=6>
      // ##FILTER=<ID=DISC,Description="Mendelian or duplicate genotype
      // discordance is high (3/5% or more)",IDX=7>
      // ...
      // ##INFO=<ID=SVM,Number=1,Type=Float,Description="Milk-SVM score for
      // variant quality, passing -0.5 or greater,IDX=6">

      if (parser.idx == (int)header_id->size()) {
        header_id->push_back(parser.id);
        header_number->push_back(parser.number);
        header_type->push_back(parser.type);
        header_desc->push_back(parser.desc);
      } else if (parser.idx < (int)header_id->size()) {
        printf("BCF index (IDX=) is reused for [%s] with IDX=%d\n",
               parser.id.c_str(), parser.idx);
        (*header_id)[parser.idx] = parser.id;
        (*header_number)[parser.idx] = parser.number;
        (*header_type)[parser.idx] = parser.type;
        (*header_desc)[parser.idx] = parser.desc;
      } else {
        printf("BCF index is invalid for [%s] with IDX=%d, skipped!\n",
               parser.id.c_str(), parser.idx);
      }
    } else {
      // do nothing
    }
  }
  printf("Total contig parse = %d, total header index used = %d\n",
         (int)header_contig_id->size(), (int)header_id->size());
  return 0;
}

int readOneInteger(BGZF* fp, int* len) {
  uint8_t val_type;
  // read the next integer
  if (1 != bgzf_read(fp, &val_type, 1)) {
    fprintf(stderr, "Wrong read!\n");
    exit(1);
  }
  switch (val_type & 0x0F) {
    case 1:  // 8 bit int
      int8_t tmp8;
      if (1 != bgzf_read(fp, &tmp8, sizeof(int8_t))) {
        *len = tmp8;
      }
      break;
    case 2:  // 16 bit int
      int16_t tmp16;
      if (1 != bgzf_read(fp, &tmp16, sizeof(int16_t))) {
        *len = tmp16;
      }
      break;
    case 3:  // 32 bit int
      int32_t tmp32;
      if (1 != bgzf_read(fp, &tmp32, sizeof(int32_t))) {
        *len = tmp32;
      }
      break;
    default:
      fprintf(stderr, "Wrong type!\n");
      exit(1);
  }
  if (val_type >> 4 != 1) {
    fprintf(stderr, "Wrong array dimension!\n");
    exit(1);
  }
  return 0;
}

// read one or two bytes of given @param type from @param fp, and report array
// $param len
int readArray(BGZF* fp, const int type, int* len) {
  uint8_t val_type;
  if (1 != bgzf_read(fp, &val_type, 1)) {
    fprintf(stderr, "Wrong read!\n");
    exit(1);
  }
  if ((val_type & 0x0F) != type) {
    fprintf(stderr, "Wrong type %d != %d!\n", val_type & 0x0F, type);
    exit(1);
  }
  uint8_t val_len = (val_type >> 4);
  if (val_len == 0) {  // missing
    *len = 0;
  } else if (val_len == 15) {  // overflowed
    readOneInteger(fp, len);
  } else {
    *len = val_len;
  }
  return 0;
}

int readString(BGZF* fp, std::string* ret) {
  int len;
  if (readArray(fp, 7, &len)) {
    fprintf(stderr, "Wrong read array!\n");
    exit(1);
  }
  printf("len of string = %d\n", len);
  ret->resize(len);
  if ((ssize_t)(len * sizeof(char)) !=
      bgzf_read(fp, (void*)ret->data(), len * sizeof(char))) {
    fprintf(stderr, "Wrong read string!\n");
    exit(1);
  }
  return 0;
}

int readInt(BGZF* fp, std::vector<int8_t>* ret) {
  int len;
  if (readArray(fp, 1, &len)) {  // 1 means 8bit integer
    fprintf(stderr, "Wrong read array!\n");
    exit(1);
  }
  printf("len of int = %d\n", len);
  ret->resize(len);
  if ((ssize_t)(len * sizeof(int8_t)) !=
      bgzf_read(fp, (void*)ret->data(), len * sizeof(int8_t))) {
    fprintf(stderr, "Wrong read string!\n");
    exit(1);
  }
  return 0;
}

int readInt(BGZF* fp, std::vector<int16_t>* ret) {
  int len;
  if (readArray(fp, 2, &len)) {  // 3 means 32bit integer
    fprintf(stderr, "Wrong read array!\n");
    exit(1);
  }
  printf("len of int = %d\n", len);
  ret->resize(len);
  if ((ssize_t)(len * sizeof(int16_t)) !=
      bgzf_read(fp, (void*)ret->data(), len * sizeof(int16_t))) {
    fprintf(stderr, "Wrong read string!\n");
    exit(1);
  }
  return 0;
}

int readFloat(BGZF* fp, std::vector<float>* ret) {
  int len;
  if (readArray(fp, 5, &len)) {  // 3 means float
    fprintf(stderr, "Wrong read array!\n");
    exit(1);
  }
  printf("len of int = %d\n", len);
  ret->resize(len);
  if ((ssize_t)(len * sizeof(float)) !=
      bgzf_read(fp, (void*)ret->data(), len * sizeof(float))) {
    fprintf(stderr, "Wrong read float!\n");
    exit(1);
  }
  return 0;
}

int readInt(BGZF* fp, std::vector<int>* ret) {
  // read one or two bytes of integers (int8_t, int16_t or int32_t) from @param
  // fp, and report array $param len
  // int readIntArray(BGZF* fp, int* len) {
  uint8_t val_type;
  if (1 != bgzf_read(fp, &val_type, 1)) {
    fprintf(stderr, "Wrong read!\n");
    exit(1);
  }
  if ((val_type & 0x0F) < 1 ||
      (val_type & 0x0F) >
          3) {  // 1, 2, 3 represents int8_t, int16_t and int32_t
    fprintf(stderr, "Wrong int type %d !\n", val_type & 0x0F);
    exit(1);
  }
  uint8_t val_len = (val_type >> 4);
  int len;
  if (val_len == 0) {  // missing
    len = 0;
  } else if (val_len == 15) {  // overflowed
    readOneInteger(fp, &len);
  } else {
    len = val_len;
  }

  ret->resize(len);
  int8_t* p8;
  int16_t* p16;
  int32_t* p32;
  switch (val_type & 0x0F) {
    case 1:
      p8 = new int8_t[len];
      if ((ssize_t)(len * sizeof(int8_t)) !=
          bgzf_read(fp, (void*)ret->data(), len * sizeof(int8_t))) {
        fprintf(stderr, "Wrong read int8_t!\n");
        exit(1);
      }
      std::copy(p8, p8 + len, ret->begin());
      delete[] p8;
      break;
    case 2:
      p16 = new int16_t[len];
      if ((ssize_t)(len * sizeof(int16_t)) !=
          bgzf_read(fp, (void*)ret->data(), len * sizeof(int16_t))) {
        fprintf(stderr, "Wrong read int16_t!\n");
        exit(1);
      }
      std::copy(p16, p16 + len, ret->begin());
      delete[] p16;
      break;
    case 3:
      p32 = new int32_t[len];
      if ((ssize_t)(len * sizeof(int32_t)) !=
          bgzf_read(fp, (void*)ret->data(), len * sizeof(int32_t))) {
        fprintf(stderr, "Wrong read int32_t!\n");
        exit(1);
      }
      std::copy(p32, p32 + len, ret->begin());
      delete[] p32;
      break;
    default:
      fprintf(stderr, "Wrong int type %d !\n", val_type & 0x0F);
      exit(1);
      break;
  }
  return 0;
}

int main(int argc, char** argv) {
  const char* fn = argv[1];
  std::string indexFile = "unnamed.scIndex";
  if (argc == 3) {
    indexFile = argv[2];
  } else {
    // indexFile = "tmp.scIndex";
  }
  printf("indexFile = %s\n", indexFile.c_str());

  BGZF* fp = bgzf_open(fn, "rb");
  if (!fp) {
    exit(1);
  }

  // read header
  char magic[5];
  if (5 != bgzf_read(fp, magic, 5)) {
    exit(1);
  }
  if (!(magic[0] == 'B' && magic[1] == 'C' && magic[2] == 'F' &&
        magic[3] == 2 && (magic[4] == 1 || magic[4] == 2))) {
    exit(1);
  }

  uint32_t l_text;
  if (4 != bgzf_read(fp, &l_text, 4)) {
    exit(1);
  }
  printf("l_text = %d\n", l_text);

  std::string s;
  int64_t bgzf_offset_before_header =
      bgzf_tell(fp);  // the beginning of header block
  s.resize(l_text);
  if (bgzf_read(fp, (void*)s.data(), l_text) != l_text) {
    fprintf(stderr, "Read failed!\n");
  }
  if (parseHeader(s, &header_contig_id, &header_id, &header_number,
                  &header_type, &header_description)) {
    fprintf(stderr, "Parse header failed!\n");
    exit(1);
  }

  int64_t bgzf_offset_after_header = bgzf_tell(fp);  // the end of header block
  size_t ptr_chrom_line = s.find("#CHROM");  // the index of "#CHROM", also the
                                             // size between beginning of header
                                             // to '#CHROM'
  if (ptr_chrom_line == std::string::npos) {
    fprintf(stderr, "Cannot find the \"#CHROM\" line!\n");
    exit(1);
  }
  printf("offset_header = %d\n", (int)ptr_chrom_line);

  bgzf_seek(fp, bgzf_offset_before_header, SEEK_SET);
  s.resize(ptr_chrom_line);
  int64_t before_chrom_size = bgzf_read(fp, (void*)s.data(), ptr_chrom_line);
  int64_t bgzf_offset_before_chrom = bgzf_tell(fp);  // the offset to #CHROM
  s.resize(l_text - before_chrom_size);
  int64_t after_chrom_size =
      bgzf_read(fp, (void*)s.data(), l_text - before_chrom_size);

  std::vector<std::string> fd;
  stringTokenize(s, "\t", &fd);
  printf("fd.size() = %d\n", (int)fd.size());
  const int64_t num_sample = (int)fd.size() - 9;
  printf("sample size = %ld\n", num_sample);
  printf("last character is s[after_chrom_size-1] = %d\n",
         s[after_chrom_size - 1]);

  if (bgzf_offset_after_header != bgzf_tell(fp)) {
    fprintf(stderr, "Messed up bgzf header\n");
    exit(1);
  }

  FILE* fIndex = fopen(indexFile.c_str(), "wb");
  int64_t num_marker = 0;
  int64_t pos = 0;
  fwrite(&num_sample, sizeof(int64_t), 1, fIndex);
  fwrite(&num_marker, sizeof(int64_t), 1, fIndex);
  fwrite(&pos, sizeof(int64_t), 1, fIndex);
  fwrite(&bgzf_offset_before_chrom, sizeof(int64_t), 1, fIndex);

  uint32_t l_shared;
  uint32_t l_indiv;
  std::vector<char> data;
  int64_t offset;

  do {
    offset = bgzf_tell(fp);
    if (4 != bgzf_read(fp, &l_shared, sizeof(uint32_t))) {
      break;  // fprintf(stderr, "Wrong read!\n"); exit(1);
    }
    if (4 != bgzf_read(fp, &l_indiv, sizeof(uint32_t))) {
      break;  // fprintf(stderr, "Wrong read!\n"); exit(1);
    }
    data.resize(l_shared + l_indiv);
    if (l_shared + l_indiv !=
        bgzf_read(fp, data.data(), (l_shared + l_indiv) * sizeof(char))) {
      break;  // fprintf(stderr, "Wrong read!\n"); exit(1);
    }
    memcpy(&pos, data.data() + 4, 4);
    fwrite(&pos, sizeof(int64_t), 1, fIndex);
    fwrite(&offset, sizeof(int64_t), 1, fIndex);

    num_marker++;
    if (num_marker % 10000 == 0) {
      printf("\rprocessed %ld markers", num_marker);
    }
  } while (true);

  if (fseek(fIndex, 0, SEEK_SET)) {
    fprintf(stderr, "fseek failed\n!");
  }
  fwrite(&num_sample, sizeof(int64_t), 1, fIndex);
  fwrite(&num_marker, sizeof(int64_t), 1, fIndex);
  fclose(fIndex);
  printf("Indexing finished with %ld samples and %ld markers\n", num_sample,
         num_marker);

  return 0;

  // parse one variant
  for (int var_idx = 0; var_idx < 2; ++var_idx) {
    printf("------------------------- var %d -------------------------\n",
           var_idx);
    uint32_t l_shared;
    if (4 != bgzf_read(fp, &l_shared, sizeof(uint32_t))) {
      fprintf(stderr, "Wrong read!\n");
      exit(1);
    }
    uint32_t l_indiv;
    if (4 != bgzf_read(fp, &l_indiv, sizeof(uint32_t))) {
      fprintf(stderr, "Wrong read!\n");
      exit(1);
    }
    printf("l_shared = %u, l_indiv = %u\n", l_shared, l_indiv);

    int32_t chrom;
    int32_t pos;  // 0-based index
    int32_t rlen;
    float qual;
    uint32_t n_allele_info;  // n_allele<<16 | n_info
    uint32_t n_fmt_sample;   // n_fmt<<24 | n_sample

    if (4 != bgzf_read(fp, &chrom, sizeof(int32_t))) {
      fprintf(stderr, "Wrong read!\n");
      exit(1);
    }
    if (4 != bgzf_read(fp, &pos, sizeof(int32_t))) {
      fprintf(stderr, "Wrong read!\n");
      exit(1);
    }
    if (4 != bgzf_read(fp, &rlen, sizeof(int32_t))) {
      fprintf(stderr, "Wrong read!\n");
      exit(1);
    }
    if (4 != bgzf_read(fp, &qual, sizeof(float))) {
      fprintf(stderr, "Wrong read!\n");
      exit(1);
    }
    if (4 != bgzf_read(fp, &n_allele_info, sizeof(uint32_t))) {
      fprintf(stderr, "Wrong read!\n");
      exit(1);
    }
    if (4 != bgzf_read(fp, &n_fmt_sample, sizeof(uint32_t))) {
      fprintf(stderr, "Wrong read!\n");
      exit(1);
    }
    int32_t n_allele = n_allele_info >> 16;
    int32_t n_info = n_allele_info & ((1 << 16) - 1);
    int32_t n_fmt = n_fmt_sample >> 24;
    int32_t n_sample = n_fmt_sample & ((1 << 24) - 1);
    printf(
        "%d:%d, rlen = %d, qual = %g, allele = %d, info = %d, fmt = %d, sample "
        "= %d\n",
        chrom, pos, rlen, qual, n_allele, n_info, n_fmt, n_sample);

    // ID: Variant identifier, typed str
    std::string var_id;
    if (readString(fp, &var_id)) {
      fprintf(stderr, "Wrong read!\n");
      exit(1);
    }
    printf("var_id = %s\n", var_id.c_str());
    // read alleles
    for (int i = 0; i < n_allele; ++i) {
      if (readString(fp, &var_id)) {
        fprintf(stderr, "Wrong read!\n");
        exit(1);
      }
      printf("allele[%d] = %s\n", i, var_id.c_str());
    }

    // FILTER
    std::vector<int8_t> filt;
    if (readInt(fp, &filt)) {
      fprintf(stderr, "Wrong read!\n");
      exit(1);
    }
    printf("filt = ");
    for (size_t i = 0; i != filt.size(); ++i) {
      printf("%d (%s) ", filt[i], header_id[filt[i]].c_str());
    }
    printf("\n");

    // info_key, info_value
    for (int i = 0; i < 5; ++i) {
      std::vector<int8_t> info_key;
      readInt(fp, &info_key);
      printf(
          "info key = %d [%s], number = [%s], type = [%s], description = [%s] "
          "\n",
          info_key[0], header_id[info_key[0]].c_str(),
          header_number[info_key[0]].c_str(), header_type[info_key[0]].c_str(),
          header_description[info_key[0]].c_str());
      std::vector<float> val_float;
      std::vector<int> val_int;

      if (header_type[info_key[0]] == "Float") {
        readFloat(fp, &val_float);
        printf("value = %g\n", val_float[0]);
      } else if (header_type[info_key[0]] == "Integer") {
        readInt(fp, &val_int);
        printf("value = %d\n", val_int[0]);
      }
    }

    std::vector<int8_t> format_key;
    readInt(fp, &format_key);
    printf(
        "format key = %d [%s], number = [%s], type = [%s], description = [%s] "
        "\n",
        format_key[0], header_id[format_key[0]].c_str(),
        header_number[format_key[0]].c_str(),
        header_type[format_key[0]].c_str(),
        header_description[format_key[0]].c_str());
    int8_t format_type;
    bgzf_read(fp, &format_type, 1);
    printf("format type = 0x%0x ", format_type);
    int format_len_per_indv =
        (format_type >> 4) *  // (num of types per indv)
        (format_type &
         ((1 << 4) -
          1));  // (bytes per type, e.g. 1 for int8_t, which is 1 byte
    printf("format len per indv = %d\n", format_len_per_indv);
    // read and discard genotypes
    int8_t* p_discard = new int8_t[format_len_per_indv];
    std::map<int16_t, int> freq;
    for (int i = 0; i < num_sample; ++i) {
      if (format_len_per_indv !=
          bgzf_read(fp, p_discard, format_len_per_indv)) {
        fprintf(stderr, "Cannot read data!\n");
        exit(1);
      }
      int16_t v = (p_discard[0] << 8) | p_discard[1];
      freq[v]++;
    }
    for (std::map<int16_t, int>::iterator iter = freq.begin();
         iter != freq.end(); ++iter) {
      if (iter->first == 0) {
        printf("./.  %d\n", iter->second);
      } else {
        printf("%d/%d %d\n", (iter->first >> 9) - 1,
               ((iter->first & ((1 << 8) - 1)) >> 1) - 1, iter->second);
      }
    }
    delete[] p_discard;
  }

#if 0
  kstring_t* str;
  str = (kstring_t*)calloc(1, sizeof(kstring_t));
  kstring_t& s = *str;
  FILE* fIndex = fopen("test.idx", "wb");
  int ret;
  int64_t offset;
  int64_t pos;
  do {
    offset = bgzf_tell(fp);
    ret = bgzf_getline(fp, '\n', &s);

    if (ret <= 0) {
      break;
    }
    size_t beg = 0;
    for (; beg < s.l; ++beg) {
      if (s.s[beg] == '\t') {
        pos = strtol(s.s + beg + 1, NULL, 0);
        break;
      }
    }
    if (s.s[0] == '#') {
      if (s.s[1] == '#') {
        continue;
      } else if (s.s[1] == '#') {  // header line
        pos = 0;
      } else {
        fprintf(stderr, "Strange header line!\n");
      }
    }
    // printf("%ld %ld\n", pos, offset);
    fwrite(&pos, sizeof(int64_t), 1, fIndex);
    fwrite(&offset, sizeof(int64_t), 1, fIndex);
  } while (1);

  bgzf_close(fp);
  fclose(fIndex);
  return 0;
#endif
}
