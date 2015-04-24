#include "ParRegion.h"

int main(int argc, char *argv[])
{
  {
    ParRegion p;
    assert(!p.isParRegion("1", 100));
    assert(!p.isParRegion("X", 60000));
    assert(p.isParRegion("X", 60001));
    assert(!p.isParRegion("X", 154931043));
    assert(p.isParRegion("X", 154931044));
    assert(p.isParRegion("23", 154931044));
    assert(!p.isParRegion("23", 155260561));
  }
  {
    ParRegion p("X", "b36");
    assert(!p.isParRegion("1", 100));
    
    assert(!p.isParRegion("X", 0));
    assert(p.isParRegion("X", 1));
    assert(p.isParRegion("X", 2709520));
    assert(!p.isParRegion("23", 2709521));
    assert(p.isParRegion("X", 154584238));
    assert(p.isParRegion("X", 154913754));
    assert(!p.isParRegion("23", 154913754));
    assert(!p.isParRegion("23", 154913755));
  }

  {
    // fake example
    ParRegion p("1", "-100");
    assert(p.getRegion().empty());
  }

  {
    // fake example
    ParRegion p("1", "100-200");
    // fprintf(stderr, "%s\n", p.getRegion().c_str());
    assert(p.getRegion() == "100-200");
  }

  {
    // fake example
    ParRegion p("1", "100-");
    assert(p.isParRegion("1", 100));
    assert(p.isParRegion("1", 59373566)); // 59373566 is chrom X length in GrCh37
  }
  return 0;
}
