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

    m.clear();
    assert (m[0][0] == 0.0);
    assert (m[2][1] == 0.0);
    
    std::vector<double> col(3, 1.0);
    m.appendCol(col);
    assert (m[0][0] == 0.0);
    assert (m[2][2] == 1.0);

    m.resize(2,3);
    std::vector<double> row(3, 2.0);
    m.appendRow(row);
    assert (m[0][0] == 0.0);
    assert (m[2][2] == 2.0);

    
    return 0;
}
