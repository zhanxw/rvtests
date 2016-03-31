#include "SimpleMatrix.h"
#include <cassert>

int main(int argc, char *argv[])
{
    SimpleMatrix m;
    assert (m.nrow() == 0);
    assert (m.ncol() == 0);

    m.resize(3, 2);
    assert (m.nrow() == 3);
    assert (m.ncol() == 2);

    m.zero();
    assert (m[0][0] == 0.0);
    assert (m[2][1] == 0.0);
    
    std::vector<double> col(3, 1.0);
    m.appendCol(col);
    assert (m[0][0] == 0.0);
    assert (m[2][2] == 1.0);

    m.resize(2,3);
    std::vector<double> row(3, 22.0);
    m.appendRow(row);
    assert (m[0][0] == 0.0);
    assert (m[2][2] == 22.0);
    m.appendRow(row);
    assert (m[3][2] == 22.0);
    
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
        
    return 0;
}
