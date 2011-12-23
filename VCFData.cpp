#include "VCFData.h"

// outputs:
// prefix + ".geno": raw genotype
// prefix + ".cgeno": raw genotype
// prefix + ".cov": raw genotype
// prefix + ".pheno": raw genotype
void VCFData::writeRawData(const char* prefix){
    std::string p = prefix;
    if (p == "") {
        p = "rvtest.raw";
    }
    this->writeGenotype( (p + ".geno").c_str());
    // this->writeCollapsedGenotype( (p + ".cgeno").c_str());
    this->writeCovariate( (p + ".cov").c_str());
    this->writePhenotype( (p + ".pheno").c_str());
};

// write genotype in R format:
//   have header: PeopleID MarkerName[0] MarkerName[1]...
//   don't have row names, e.g. 1, 2, 3, ...
//   genotype are people x marker
void VCFData::writeGenotype(const char* fn){
    this->writeTable(fn, this->genotype, this->marker2Idx, this->people2Idx, "MarkerName");
};
// void VCFData::writeCollapsedGenotype( const char* fn){
//     this->writeTable(fn, this->collapsedGenotype, this->people2Idx, this->set2Idx, "PeopleID");
// };
/** write covariate to file, format is as following:
 *  header line: PeopleName CovName1 CovName2 ...
 *  content line: P1 1.0 2.0 ...
 */
void VCFData::writeCovariate(const char* fn) {
    this->writeTable(fn, this->covariate, this->people2Idx, this->covariate2Idx, "PeopleID");
};
void VCFData::writePhenotype(const char* fn) {
    this->writeTable(fn, this->phenotype, this->people2Idx, this->phenotype2Idx, "PeopleID");
};

/**
 * write @param data in to a R-readable format
 */
void VCFData::writeTable(const char* fn,
                         Matrix* data, OrderedMap<std::string, int> & rowName,
                         OrderedMap<std::string, int> & colName,
                         const char* upperLeftName) {
    if (!data || data->rows == 0 || data->cols == 0)
        return;

    if (rowName.size() != data->rows) {
        fprintf(stderr, "Row number does not match!");
        return;
    }
    if (colName.size() != data->cols) {
        fprintf(stderr, "Col number does not match!");
        return;
    }

    FileWriter fw(fn);
    // header
    fw.write(upperLeftName);
    for (int i = 0; i < colName.size(); i++){
        fw.write("\t");
        if (colName.size() == 0){
            fw.write("\".\"");
        }else {
            fw.write("\"");
            fw.write(colName.keyAt(i).c_str());
            fw.write("\"");
        }
    };
    fw.write("\n");

    // content
    int numCol = colName.size();
    int numRow = rowName.size();
    for (int r = 0; r < numRow; r++ ){
        fw.write(rowName.keyAt(r).c_str());
        for (int c = 0; c < numCol; c++) {
            fw.printf("\t%.3f", (double)((*data)[r][c]));
        }
        fw.write("\n");
    }
    fw.close();
};

/**
 * read @param fn into @param data from R-readable format
 * REQURIE:
 *   @param rowName and @param colName are both treated like string
 *   @param data should be integer/float number
 *   @param defaultValue, when data cannot be converted to double, use this value instead.
 */
int VCFData::readTable(const char* fn,
                         Matrix* data, OrderedMap<std::string, int> * rowName,
                         OrderedMap<std::string, int> * colName,
                        std::string* upperLeftName,
                        double defaultValue) {
    if (!data || data->rows == 0 || data->cols == 0 || !upperLeftName)
        return -1;

    int invalidConversion = 0;
    LineReader lr(fn);
    std::vector <std::string> fd;
    std::vector <std::string> header;
    int lineNo = 0;
    int nCol = -1; // col number including row name.
    while (lr.readLineBySep(&fd, " \t")){
        if (lineNo == 0){
            header = fd;
            // since the first column may or may not have header, 
            // we will later check the size of header and set the values of colName
            continue;
        }
        if (lineNo == 1){
            nCol = fd.size();
            if (fd.size() == header.size()) {
                // first column have header
                for (int r = 1; r < header.size(); r++) 
                    (*rowName)[header[r]] = r - 1;
            } else{
                // first column does not have header
                for (int r = 0; r < header.size(); r++) 
                    (*rowName)[header[r]] = r ;
            }
        }

        if (fd.size() != nCol) {
            fprintf(stderr, "Inconsistent column number at line %d, skipping...\n", lineNo);
            continue;
        }

        int row = lineNo - 1;
        data->Dimension(lineNo, nCol - 1);
        (*rowName)[fd[0]] = row;
        for (int col = 1; col < nCol; col++){
            if (!str2double( fd[col].c_str(),  &(*data)[row][col])){
                (*data)[row][col] = defaultValue;
                invalidConversion++;
            }
        }
        lineNo++;
    };
    return invalidConversion;
};

int VCFData::readPlinkTable(const char* fn,
                            Matrix* data, 
                            OrderedMap<std::string, int> * rowName,
                            OrderedMap<std::string, int> * colName,
                            double defaultValue){
    if (!fn || !data)
        return -1;

    assert(rowName && colName);
    rowName->clear();
    colName->clear();

    int invalidConversion = 0;
    LineReader lr(fn);
    std::vector <std::string> fd;
    int lineNo = 0;
    int nCol = -1; // col number including row name.
    while (lr.readLineBySep(&fd, " \t")){
        if (lineNo == 0){
            nCol = fd.size();
            for (int c = 2; c < fd.size(); c++) {// skip first 2 column
                (*colName)[fd[c]] = c - 1;
            }
            continue;
        }
        if (fd.size() != nCol) {
            fprintf(stderr, "Inconsistent column number at line %d, skipping...\n", lineNo);
            continue;
        }
        int row = rowName->size(); 

        if (rowName->find(fd[1])){
            fprintf(stderr, "Duplicate sample: %s\n", fd[1].c_str());
        } else {
            data->Dimension(row + 1, nCol - 1);
        }

        (*rowName)[fd[1]] = row;
        for (int col = 2; col < nCol; col++){
            if (!str2double( fd[col].c_str(),  &(*data)[row][col])){
                (*data)[row][col] = defaultValue;
                invalidConversion++;
            }
        }
        lineNo++;
    };
    return invalidConversion;
};

/**
 * @return -1: error
 *          >=0 : num of individuals skipped
 */
int VCFData::readPlinkPhenotypeSkipMissing(const char* fn, const char* selectedCol,
                                           Matrix* data, 
                                           OrderedMap<std::string, int> * rowName,
                                           OrderedMap<std::string, int> * colName){
    if (!fn || !data)
        return -1;

    if (strlen(fn) == 0) {
        fprintf(stderr, "Emptye phenotpye file!\n");
        return -1;
    }

    assert(rowName && colName);
    rowName->clear();
    colName->clear();

    int skipped = 0;
    LineReader lr(fn);
    std::vector <std::string> fd;
    int lineNo = 0;
    int nCol = -1; // col number including row name.
    while (lr.readLineBySep(&fd, " \t")){
        if (lineNo == 0){
            nCol = fd.size();
            if (nCol < 2) {
                fprintf(stderr, "Too few column for pheontype file: %s.\n", fn);
                return -1;
            }
            if (selectedCol  && strlen(selectedCol) > 0) {
                if (fd[0] == "FID" && fd[1] == "IID") {
                    for (int r = 2; r < fd.size(); r++) {// skip first 2 column
                        (*colName)[fd[r]] = r - 1;
                    }
                    continue;
                } else{
                    fprintf(stderr, "Phenotype file header should be FID and IID: %s.\n", fn);                    
                    return -1;
                }
            } else {
                // read only 3rd column
                // TODO: "add support to selected column"
                // for (int r = 2; r < fd.size(); r++) {// skip first 2 column
                //     (*rowName)[fd[r]] = r - 1;
                // }
                (*colName)["Pheno"] = 0;
            }
        }
        if (fd.size() != nCol) {
            fprintf(stderr, "Inconsistent column number at line %d, skipping...\n", lineNo);
            continue;
        }
        int row = rowName->size() ;

        if (rowName->find(fd[1])){
            fprintf(stderr, "Duplicate sample: %s\n", fd[1].c_str());
        } else{
            data->Dimension(row + 1, 1);
        }

        int col = 2; // use 3rd column as phenotype
        if (selectedCol && colName->find(selectedCol))
            col = (*colName)[selectedCol];
        if (!str2double( fd[col].c_str(),  &(*data)[row][0])){
            fprintf(stderr, "Skiping missing phenotype of individual %s.\n", fd[1].c_str());
            skipped ++;
            continue;
        } else {
            (*rowName)[fd[1]] = row;
        };
        lineNo ++;
    };
    return skipped;
};


