#include <cassert>
#include <cmath>
#include "SimpleMatrix.h"
int main(int argc, char* argv[]) {
  {
    SimpleMatrix m;
    assert(m.nrow() == 0);
    assert(m.ncol() == 0);

    m.resize(3, 2);
    assert(m.nrow() == 3);
    assert(m.ncol() == 2);

    m.zero();
    assert(m[0][0] == 0.0);
    assert(m[2][1] == 0.0);

    std::vector<double> col(3, 1.0);
    m.appendCol(col);
    assert(m[0][0] == 0.0);
    assert(m[2][2] == 1.0);

    m.resize(2, 3);
    std::vector<double> row(3, 22.0);
    m.appendRow(row);
    assert(m[0][0] == 0.0);
    assert(m[2][2] == 22.0);
    m.appendRow(row);
    assert(m[3][2] == 22.0);

    const char* fn = "tmp.mat.out";
    assert(m.writeFile(fn) == 0);
    SimpleMatrix m2;
    assert(m2.readFile(fn) == 0);
    assert(m2.nrow() == m.nrow());
    assert(m2.ncol() == m.ncol());
    assert(m2[0][0] == m[0][0]);
    assert(m2[1][2] == m[1][2]);
    assert(m2[3][2] == m[3][2]);

    // m is a 4 by 3 matrix
    m.deleteRow(0);
    m.deleteCol(0);
    m.deleteRow(m.nrow() - 1);
    m.deleteCol(m.ncol() - 1);
    assert(m.nrow() == 2);
    assert(m.ncol() == 1);
  }

  {
    // m:
    // 0 1 2
    // 3 4 5
    SimpleMatrix m;
    m.resize(2, 3);
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 3; ++j) {
        m[i][j] = i * 3 + j;
      }
    }

    std::vector<int> ind = {-1, 0, 1, 1};
    m.reorderRow(ind);
    // now m should be
    // m:
    // 0 1 2
    // 3 4 5
    // 3 4 5
    assert(m.getRowName()[0] == "R1");
    assert(m.getRowName()[2] == "R2");
    assert(m.nrow() == 3 && m.ncol() == 3);
    assert(m[2][0] == 3);
    assert(m[2][2] == 5);
  }

  {
    // m:
    // 0 1 2
    // 3 4 5
    SimpleMatrix m;
    m.resize(2, 3);
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 3; ++j) {
        m[i][j] = i * 3 + j;
      }
    }

    std::vector<std::string> ind = {"XX", "R1", "R2", "R2"};
    m.reorderRow(ind);
    // now m should be
    // m:
    // 0 1 2
    // 3 4 5
    // 3 4 5
    assert(m.getRowName()[0] == "R1");
    assert(m.getRowName()[2] == "R2");
    assert(m.nrow() == 3 && m.ncol() == 3);
    assert(m[2][0] == 3);
    assert(m[2][2] == 5);
  }

  {
    // m:
    // 0 1 2
    // 3 4 5
    SimpleMatrix m;
    m.resize(2, 3);
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 3; ++j) {
        m[i][j] = i * 3 + j;
      }
    }

    std::vector<std::string> ind = {"New1", "New2", "New3"};
    m.addRow(ind, NAN);
    assert(std::isnan(m[2][0]));
    assert(std::isnan(m[2][2]));
    assert(std::isnan(m[4][0]));
    assert(std::isnan(m[4][2]));
    assert(m.getRowName()[2] == "New1");
    assert(m.getRowName()[4] == "New3");
    assert(m.nrow() == 5 && m.ncol() == 3);

    // now m should be
    // m:
    // 0   1   2
    // 3   4   5
    // 1.5 2.5 3.5
    // 1.5 2.5 3.5
    // 1.5 2.5 3.5
    m.imputeMissingToMeanByCol();
    assert(m[2][0] = 1.5);
    assert(m[2][2] = 3.5);
    assert(m[4][0] = 1.5);
    assert(m[4][2] = 3.5);
  }

  return 0;
}
