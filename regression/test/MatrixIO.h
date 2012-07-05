#ifndef _MATRIXIO_H_
#define _MATRIXIO_H_

void LoadVector(const char* fn, Vector& v){
    LineReader lr(fn);
    std::vector< std::string> s;
    int lineNo = 0;
    while (lr.readLineBySep(&s, " \t")) {
        lineNo ++;
        v.Dimension(lineNo);
        v[lineNo - 1] = atof(s[0].c_str());
    }
};
void LoadMatrix(const char* fn, Matrix& m){
    LineReader lr(fn);
    std::vector< std::string> s;
    int lineNo = 0;
    while (lr.readLineBySep(&s, " \t")) {
        lineNo ++;
        m.Dimension(lineNo , s.size());
        for (int j = 0; j < s.size() ; j++ ) {
            m[lineNo - 1][j] = atof(s[j].c_str());
        }
    }


};

void Print(Vector& v){
    for (int i = 0; i < v.Length() ; i++){
        if (i) { fprintf(stdout, "\t");}
        fprintf(stdout, "%.3f", v[i]);
    }
};
void Print(Matrix& m){
    for (int i = 0; i < m.rows ; i++){
        if (i) { fprintf(stdout, "\t");}
        for (int j = 0; j < m.cols; j++){
            if (j) { fprintf(stdout, "\t");}
            fprintf(stdout, "%.3f", m[i][j]);
        }
    }
}
void Print(double& d){
    fprintf(stdout, "%.3f", d);
}


#endif /* _MATRIXIO_H_ */
