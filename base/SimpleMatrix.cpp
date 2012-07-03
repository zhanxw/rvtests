#include "SimpleMatrix.h"
#include "IO.h"
#include "TypeConversion.h"

/**
 * @return 0: success
 */
int SimpleMatrix::readFile(const char* f){
    LineReader lr(f);
    std::vector<std::string> fd;
    std::vector< double > d;
    while(lr.readLineBySep(&fd, " \t")){
        d.resize(fd.size());
        for (unsigned int i = 0; i < fd.size(); i++)
            d[i] = atof(fd[i]);
        if (mat.size() && d.size() != mat[mat.size() -1 ].size()) {
            fprintf(stderr, "Column width does not fit!\n");
            return -1;
        }
        this->mat.push_back(d);
    };
    return 0;
};

int SimpleMatrix::writeFile(const char* f){
    FileWriter fw(f);
    for (unsigned int i = 0; i < mat.size(); i++){
        for (unsigned int j = 0; j < mat.size(); j++) {
            fw.printf("%f", mat[i][j]);
            if (j) 
                fw.write("\t");
        }
        fw.write("\n");
    }
    return 0;
};
