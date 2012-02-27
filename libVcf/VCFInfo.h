#ifndef _VCFINFO_H_
#define _VCFINFO_H_

class VCFInfoValue{
  public:
    int fingerMark; 
    VCFValue* value;
  public:
    VCFInfoValue(){
        this->value = new VCFValue;
        assert(this->value);
    }
    ~VCFInfoValue(){
        delete this->value;
        this->value = NULL;
    }
};

class VCFInfo{
public:
    VCFInfo() {
        this->hasParsed = false;
    };
    ~VCFInfo(){
        for (this->tableIter = this->table.begin(); 
             this->tableIter != this->table.end(); 
             this->tableIter ++){
            if (this->tableIter->second != NULL){
                delete this->tableIter->second;
            }
        }
    };
    void reset() { this-> hasParsed = false;};
    const char* getTag(const char* tag) {
        if (!tag || tag[0] == '\0') 
            return NULL;

        std::string s = tag;
        if (!hasParsed){
            this->parseActual();
        }
        this->tableIter = this->table.find(s);
        if (this->tableIter != this->table.end()){
            return this->tableIter->second->value->toStr();
        } else {
            return NULL;
        }
    } ;
    void parse(VCFValue* v) {
        this->self = v;
        this->hasParsed = false;
        this->fingerMark ++ ;
    };
    void parseActual(){
        this->hasParsed = true;
        this->origString = this->self->line + this->self->beg;
        const char* line = this->self->line;
        int b = this->self->beg;
        int e = this->self->beg;
        static std::string key;

        while ( e < this->self->end){
            key.clear();
            // find tag name;
            while(line[e] != '='){
                key.push_back(this->self->line[e++]);
            }

            b = e + 1; // skip '='
        
            VCFInfoValue* f;
            this->tableIter = this->table.find(key);
            if ( this->tableIter == this->table.end()){
                f = new VCFInfoValue;
                this->table[key] = f;
                f->fingerMark = this->fingerMark;
            } else {                
                f = this->tableIter->second;
                if (f->fingerMark == this->fingerMark) {
                    fprintf(stdout, "Duplicated key [%s] in INFO field!\n", key.c_str());
                } else {
                    f->fingerMark = this->fingerMark;
                }
            }

            this->self->parseTill(';', b, f->value);
            f->value->line[f->value->end] = '\0';
            e = f->value->end + 1;
        }
    };
  public:
    const char* getOrigString() {
        if (this->hasParsed)
            return this->origString.c_str();
        else 
            return self->toStr();
    }
private:
    VCFValue* self;
    bool hasParsed;
    int fingerMark;
    std::map<std::string, VCFInfoValue*> table;
    std::map<std::string, VCFInfoValue*>::iterator tableIter;
    std::string origString;
}; // VCFInfo

#endif /* _VCFINFO_H_ */
