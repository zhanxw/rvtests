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
    this->writeCollapsedGenotype( (p + ".cgeno").c_str());
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
void VCFData::writeCollapsedGenotype( const char* fn){
    this->writeTable(fn, this->collapsedGenotype, this->people2Idx, this->set2Idx, "PeopleID");
};
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
            fw.printf("\t%d", (int)((*data)[r][c]));
        }
        fw.write("\n");
    }
    fw.close();
};

