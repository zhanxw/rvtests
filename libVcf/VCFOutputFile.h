#ifndef _VCFOUTPUTFILE_H_
#define _VCFOUTPUTFILE_H_

#include "VCFUtil.h"

class VCFOutputFile{
public:
    VCFOutputFile(const char* fn) {
        this->fp = new FileWriter(fn);
        if (!this->fp){
            REPORT("Cannot create VCF file!");
            abort();
        }
    };
    ~VCFOutputFile(){
        if (this->fp){
            delete this->fp;
            this->fp = NULL;
        }
    };
    void writeHeader(const VCFHeader* h){
        for (int i = 0; i < h->size(); i++){
            this->fp->writeLine(h->at(i).c_str());
        }
    };
    void writeRecord(VCFRecord* r){
        this->fp->printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
                         r->getChrom(),
                         r->getPos(),
                         r->getID(),
                         r->getRef(),
                         r->getAlt(),
                         r->getQual(),
                         r->getInfo(),
                         r->getFilt(),
                         r->getFormat());
        VCFPeople& p = r->getPeople();
        for (int i = 0; i < p.size() ; i ++ ) {
            VCFIndividual* indv = p[i];
            this->fp->printf("\t%s", indv->getData().toStr());
        }
        this->fp->printf("\n");
    };
private:
    FileWriter* fp;
};

#endif /* _VCFOUTPUTFILE_H_ */
