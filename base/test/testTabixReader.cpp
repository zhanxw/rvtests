#include "TabixReader.h"

int main(int argc, char *argv[])
{
  {
    const char* fn = "test.tbl.gz";
    TabixReader tr(fn);
    tr.addRange("1:2-3");
    std::string line;
    while(tr.readLine(&line)) {
      printf("%s\n", line.c_str());
    }
  }
  {
    const char* fn = "/net/fantasia/home/zhanxw/mycode/vcf2geno/metaTest/tmp2.MetaScore.assoc.gz";
    TabixReader tr(fn);
    tr.addRange("1:564765-871276");
    std::string line;
    while(tr.readLine(&line)) {
      printf("%s\n", line.c_str());
    }
  }

  return 0;
}
