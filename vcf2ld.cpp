#include "Argument.h"
#include "IO.h"
#include "tabix.h"

#include <cassert>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "Utils.h"
#include "VCFUtil.h"

// #include "MathVector.h"
// #include "MathMatrix.h"

int calculateLD(std::vector<int>& lastGeno, std::vector<int>& geno, int* n, double* d) {
    assert(lastGeno.size() == geno.size());
    assert( n && d);

    n[0] = n[1] = n[2] = n[3] = 0;
    int l = lastGeno.size();
    for (int i = 0; i < l ; i ++ ){
        assert ( 0 == lastGeno[i] || 1== lastGeno[i]);
        assert ( 0 == geno[i] || 1== geno[i]);
        int idx = lastGeno[i] * 2 + geno[i];
        n[idx] ++ ;
    }
    
    int nTotal = n[0] + n[1] + n[2] + n[3];
    if (nTotal == 0) {
        return -1;
    };
    // Aa x Bb table (last x current)
    //          0 B ,   1 b
    //   0 A     0       1
    //   1 a     2       3
    double pAB = 1.0 * n[0] / nTotal;
    double pAb = 1.0 * n[1] / nTotal;
    double paB = 1.0 * n[2] / nTotal;
    double pab = 1.0 * n[3] / nTotal;
    
    double pA = pAB + pAb;
    double pB = pAB + paB;
    *d = pAB - pA * pB;

    return 0;
};

int main(int argc, char** argv){
    time_t currentTime = time(0);
    fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

    ////////////////////////////////////////////////
    BEGIN_PARAMETER_LIST(pl)
        ADD_PARAMETER_GROUP(pl, "Input/Output")
        ADD_STRING_PARAMETER(pl, inVcf, "--inVcf", "input VCF File")
        ADD_STRING_PARAMETER(pl, outLD, "--outLD", "output file prefix")
        ADD_PARAMETER_GROUP(pl, "People Filter")
        ADD_STRING_PARAMETER(pl, peopleIncludeID, "--peopleIncludeID", "give IDs of people that will be included in study")
        ADD_STRING_PARAMETER(pl, peopleIncludeFile, "--peopleIncludeFile", "from given file, set IDs of people that will be included in study")
        ADD_STRING_PARAMETER(pl, peopleExcludeID, "--peopleExcludeID", "give IDs of people that will be included in study")
        ADD_STRING_PARAMETER(pl, peopleExcludeFile, "--peopleExcludeFile", "from given file, set IDs of people that will be included in study")
        ADD_PARAMETER_GROUP(pl, "Site Filter")
        ADD_STRING_PARAMETER(pl, rangeList, "--rangeList", "Specify some ranges to use, please use chr:begin-end format.")
        ADD_STRING_PARAMETER(pl, rangeFile, "--rangeFile", "Specify the file containing ranges, please use chr:begin-end format.")
        END_PARAMETER_LIST(pl)
        ;    

    pl.Read(argc, argv);
    pl.Status();
    
    if (FLAG_REMAIN_ARG.size() > 0){
        fprintf(stderr, "Unparsed arguments: ");
        for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++){
            fprintf(stderr, " %s", FLAG_REMAIN_ARG[i].c_str());
        }
        fprintf(stderr, "\n");
        abort();
    }

    REQUIRE_STRING_PARAMETER(FLAG_inVcf, "Please provide input file using: --inVcf");

    const char* fn = FLAG_inVcf.c_str(); 
    VCFInputFile vin(fn);

    // set range filters here
    // e.g.     
    // vin.setRangeList("1:69500-69600");
    vin.setRangeList(FLAG_rangeList.c_str());
    vin.setRangeList(FLAG_rangeFile.c_str());

    // set people filters here
    if (FLAG_peopleIncludeID.size() || FLAG_peopleIncludeFile.size()) {
        vin.excludeAllPeople();
        vin.includePeople(FLAG_peopleIncludeID.c_str());
        vin.includePeopleFromFile(FLAG_peopleIncludeFile.c_str());
    }
    vin.excludePeople(FLAG_peopleExcludeID.c_str());
    vin.excludePeopleFromFile(FLAG_peopleExcludeFile.c_str());
    
    // let's write it out.
    FILE* fout = fopen( (FLAG_outLD + ".ld").c_str(), "wt");
    assert(fout);

    std::string lastChrom, chrom;
    int lastPos, pos;
    std::string lastRs, rs;
        
    fputs("ChrA\tPosA\tMarkerA\tChrB\tPosB\tMarkerB\tN00\tN01\tN10\tN11\tDprime\n", fout);

    std::vector<int> lastGeno;
    std::vector<int> geno;
    int lineNo = 0;
    while (vin.readRecord()){
        lineNo ++;
        VCFRecord& r = vin.getVCFRecord(); 
        VCFPeople& people = r.getPeople();
        VCFIndividual* indv;

        if (lineNo == 1) {
            chrom = r.getChrom();
            pos = r.getPos();
            rs = r.getID();

            for (int i = 0; i < people.size(); i++) {
                indv = people[i];
                // assume GTidx = 0;
                const int GTidx = 0;
                int g1 = (*indv)[GTidx].getAllele1();
                int g2 = (*indv)[GTidx].getAllele2();
                geno.push_back(g1);
                geno.push_back(g2);
            }
            continue;
        }

        // swap 
        lastGeno = geno;
        lastChrom = chrom;
        lastPos = pos;
        lastRs = rs;

        // update 
        chrom = r.getChrom();
        pos = r.getPos();
        rs = r.getID();
        
        for (int i = 0; i < people.size(); i++) {
            indv = people[i];
            // assume GTidx = 0;
            const int GTidx = 0;
            int g1 = (*indv)[GTidx].getAllele1();
            int g2 = (*indv)[GTidx].getAllele2();
            geno[i * 2 ] = (g1);
            geno[i * 2 + 1] = (g2);
        }

        // output
        fprintf(fout, "%s\t%d\t%s\t", lastChrom.c_str(), lastPos, lastRs.c_str());
        fprintf(fout, "%s\t%d\t%s\t", chrom.c_str(), pos, rs.c_str());

        int n[4] = {0};
        double d = 0.0; // D prime
        calculateLD(lastGeno, geno, n, &d);
        fprintf(fout, "%d\t%d\t%d\t%d\t%.2f\n", n[0], n[1], n[2], n[3], d);
    };

    return 0; 
};
