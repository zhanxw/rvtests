#ifndef _VCFVALUE_H_
#define _VCFVALUE_H_

/**
 * the versatile format to store value for VCF file.
 */
class VCFValue{
public:
    int beg; // inclusive
    int end; // exclusive, and beg <= end
    char* line;
public:
    int toInt() const{ 
        return atoi(line+beg);
    };
    void toInt(int* i) const {
        *i = atoi(line+beg);
    };
    double toDouble() const {
        return atof(line+beg);        
    };
    void toDouble(double* d) const {
        *d = strtod(line+beg, 0);
    };
    const char* toStr() {
        return (line + beg);
        /* this->retStr.clear(); */
        /* for (int i = beg; i < end; i++){ */
        /*     this->retStr.push_back(line[i]); */
        /* } */
        /* return (this->retStr.c_str()); */
    };
    void toStr(std::string* s) const { 
        s->clear();
        for (int i = beg; i < end; i++){
            s->push_back(line[i]);
        }
    };
public:
    // try to convert to genotype
    int getGenotype(){
        int g = 0;
        int p = beg;
        if (line[p] == '.') 
            return MISSING_GENOTYPE;
        if (line[p] < '0') 
            REPORT("Wrong genotype detected. [1]");
        else 
            g += line[p] - '0';
        
        p ++ ;
        if (p == end) 
            return g;
        
        p ++;
        if (p == end)
            REPORT("Wrong genotype length = 2");
        if (line[p] == '.')
            return MISSING_GENOTYPE;
        if (line[p] < '0')
            REPORT("Wrong genotype detected. [2]");
        else 
            g += line[p] - '0';
        
        return g;
    };
    int getAllele1(){
        int g = 0;
        int p = beg;
        if (line[p] == '.') 
            return MISSING_GENOTYPE;
        if (line[p] < '0') 
            REPORT("Wrong genotype detected. [1]");
        else 
            g += line[p] - '0';
        return g;
    };
    int getAllele2(){
        int g = 0;
        int p = beg + 2;
        if (p >= end) 
            return MISSING_GENOTYPE; 
        if (line[p] == '.') 
            return MISSING_GENOTYPE;
        if (line[p] < '0') 
            REPORT("Wrong genotype detected. [2]");
        else 
            g += line[p] - '0';
        return g;
    };
    bool isPhased(){
        if (end - beg != 3) return false;
        int p = beg + 1;
        if (line[p] == '/') return false;
        if (line[p] == '|') return true;
        return false;
    };
    bool isHaploid(){
        return (end - beg == 1);
    };
    
    /**
     * @return 0, a sub-range has successfully stored in @param result
     */
    inline int parseTill(const char c, int beg, VCFValue* result){
        assert(result);
        if (beg < this->beg || beg >= this->end) return -1;

        result->line = this->line;
        result->beg = beg;
        result->end = beg;
        while( this->line[result->end] != c && result->end < this->end) {
            result->end++;
        }
        return 0;
    };
    int parseTill(const char c, VCFValue* result){
        return parseTill(c, this->beg, result);
    };
  /* private: */
  /*   std::string retStr; // this held for the return value; */
};

#endif /* _VCFVALUE_H_ */
