#include "TypeConversion.h"
#include <climits>
#include <string>

int main(int argc, char *argv[])
{
  {
    // test toStringWithComma
    int x[] = {INT_MIN, -1234567890, -123456, -12345, -1000, -999, -1,
               0, 1, 999, 1000, 12345, 123456, 1234567890, INT_MAX};
    int n = sizeof(x)/ sizeof(x[0]);
    for (int i = 0; i < n; ++i) {
      printf("%-15d : %s\n", x[i], toStringWithComma(x[i]).c_str());
    }
  }

  {
    // test toString
    std::vector<double> a;
    a.push_back(1.1);
    a.push_back(2.0);
    std::string s = floatToString(a);
    printf("%s\n", s.c_str());
  }
  return 0;
}
