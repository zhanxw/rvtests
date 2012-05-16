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
#include "StringHash.h"
#include "Error.h"
#include <deque>

#define DEBUG_Youna 0

// #include "MathVector.h"
// #include "MathMatrix.h"

#define min(x, y) ( (x) < (y) ? (x) : (y) )

int calculateLD(std::vector<int>& lastGeno, std::vector<int>& geno, int* n, double* d) {
    assert(lastGeno.size() == geno.size());
    assert( n && d);

    n[0] = n[1] = n[2] = n[3] = 0;
    int l = lastGeno.size();
    for (int i = 0; i < l ; i ++ ){
        assert ( 0 == lastGeno[i] || 1== lastGeno[i]);
        assert ( 0 == geno[i] || 1 == geno[i]);
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
    double pa = paB + pab;
    double pb = pAb + pab;

	if (DEBUG_Youna == 1) {
		printf("pAB = %1.3e, pA = %1.3e, pB = %1.3e\n", pAB, pA, pB);
	}

    // D
    d[0] = pAB - pA * pB;
    // D'  (D prime)
    double Dmax = d[0] < 0.0 ?  (  min(pA* pB, pa * pb) ) : ( min(pA * pb, pa * pB));
    d[1] = d[0] / (Dmax + 1e-30);
    d[2] = d[0] / sqrt(pA * pa * pB * pb) ;

	// treating special cases -- only one cell has nonzero count
	if (d[0] == 0){
		for (int i = 0; i < 3; i ++) {
			d[i] = 1.0;
		}
	}

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
        ADD_PARAMETER_GROUP(pl, "Range Filter")
        ADD_STRING_PARAMETER(pl, rangeList, "--rangeList", "Specify some ranges to use, please use chr:begin-end format.")
        ADD_STRING_PARAMETER(pl, rangeFile, "--rangeFile", "Specify the file containing ranges, please use chr:begin-end format.")
        ADD_PARAMETER_GROUP(pl, "Site Filter")
        ADD_STRING_PARAMETER(pl, siteID, "--siteID", "Specify the sites to be extracted from the vcf file, separated by common")
        ADD_STRING_PARAMETER(pl, siteFile, "--siteFile", "Specify the file to contain the site to be extract from the vcf file.")
        ADD_PARAMETER_GROUP(pl, "Window Size")
        ADD_INT_PARAMETER(pl, wdSize, "--wdSize","Specify the window size for the LD calculation")
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
    vin.setRangeFile(FLAG_rangeFile.c_str());

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

	fputs("ChrA\tPosA\tMarkerA\tChrB\tPosB\tMarkerB\tN00\tN01\tN10\tN11\tD\tDprime\tr\n", fout);

	std::vector<int> lastGeno;
	std::vector<int> geno;

	// create a hash to store all the sites to be included
	if (DEBUG_Youna == 1) {
		printf("siteFile = %s\n",FLAG_siteFile.c_str() );
	}
	LineReader lr(FLAG_siteFile.c_str());
	std::vector<std::string> fd;
	StringIntHash includeSiteHash;

	while (lr.readLineBySep( & fd, "\t ")){
		includeSiteHash.Add((fd[0]+":"+fd[1]).c_str(),0);
	}

	if (DEBUG_Youna == 1) {
		printf("The size of includeSiteHash is %d \n",includeSiteHash.Entries());
	}

	String siteChrPos;
	int wdSize = FLAG_wdSize;

	std::deque< std::vector<int> > geno_deque;
	std::deque<int> pos_deque; // do not need a chromosome number deque here because we do not calcualte LD for different chromosomes here

	// read the first variant
	bool found1stInList = false;
	while (vin.readRecord() && !found1stInList) {
		VCFRecord& r = vin.getVCFRecord();
		VCFPeople& people = r.getPeople();
		VCFIndividual* indv;

		// get the first
		lastChrom = r.getChrom();
		lastPos = r.getPos();
		lastRs = r.getID();

		// get the variantID in chrom:Pos format
		siteChrPos.Clear();
		siteChrPos = lastChrom.c_str();
		siteChrPos += ":";
		siteChrPos += lastPos;

		if (includeSiteHash.Find(siteChrPos) == -1) {
			continue;
		} else {
			found1stInList = true;
			// get the genotypes
			geno.resize(2*people.size());
			for (int i = 0; i < people.size(); i++) {
				indv = people[i];
				const int GTidx = 0;
				int g1 = (*indv)[GTidx].getAllele1();
				int g2 = (*indv)[GTidx].getAllele2();
				geno[i * 2] = g1;
				geno[i *2 + 1] = g2;
			}

			geno_deque.push_back(geno);
			pos_deque.push_back(lastPos);
			break;
		}
	}

	if (!found1stInList) { // if still cannot find the 1st in List
		error("The vcf doesn't have any variants included in the list!\n");
	}

	// starting with the second variants, read through the file
	while (vin.readRecord()) {
		siteChrPos.Clear();
		// get the first line
		VCFRecord& r = vin.getVCFRecord();
		VCFPeople& people = r.getPeople();
		VCFIndividual* indv;

		// get the chrom and position number
		chrom = r.getChrom();
		pos = r.getPos();
		rs = r.getID();

		// get the variantID in chrom:Pos format
		siteChrPos.Clear();
		siteChrPos = chrom.c_str();
		siteChrPos += ":";
		siteChrPos += pos;

		if (includeSiteHash.Find(siteChrPos) == -1) {
			continue;
		}

		// get the genotype
		for (int i = 0; i < people.size(); i++) {
			indv = people[i];
			const int GTidx = 0;
			int g1 = (*indv)[GTidx].getAllele1();
			int g2 = (*indv)[GTidx].getAllele2();
			geno[i * 2] = g1;
			geno[i * 2 + 1] = g2;
		}

		if (DEBUG_Youna == 1) {
			for (int i = 0; i < pos_deque.size(); i ++) {
				printf("pos_deque[%d] = %d  ", i, pos_deque[i]);
			}
			printf("\n");
		}

		// read in the genotypes
		if ( (pos - lastPos) <= wdSize) { // if the current variant is in scope, add the genotypes in the deque
			geno_deque.push_back(geno);
			pos_deque.push_back(pos);
		} else { // if the current variant is not in scope,
			// continue to calcualte the LD in the existing deque until find the last postion that falls in scope with the current variant
			while (geno_deque.size() > 0 && (pos - lastPos) > wdSize) { // calcualte LD until the first in scope variant
				if (geno_deque.size() > 1) { // if there are at least two variants in the deque, we can do LD calculation
					for (int j = 1; j < geno_deque.size(); j ++) {
						fprintf(fout, "%s\t%d\t", lastChrom.c_str(), pos_deque[0]);
						fprintf(fout, "%s\t%d\t", chrom.c_str(), pos_deque[j]);

						int n[4] = {0};
						double d[3] = {0.0}; // D prime
						if (DEBUG_Youna == 1) {
							printf("pos_deque[%d] = %d, pos_deque[%d] = %d\n", 0, pos_deque[0], j, pos_deque[j]);
						}
						calculateLD(geno_deque[0], geno_deque[j], n, d);
						fprintf(fout, "%d\t%d\t%d\t%d\t", n[0], n[1], n[2], n[3]);
						fprintf(fout, "%.6lf\t%.6lf\t%.6lf\n", d[0], d[1], d[2]);
					}
				}
				geno_deque.pop_front();
				pos_deque.pop_front();
				if (geno_deque.size() > 0) {
					lastPos = pos_deque[0];
				}
			}

			geno_deque.push_back(geno);
			pos_deque.push_back(pos);
			lastPos = pos_deque[0];
		}
	}

	if (DEBUG_Youna == 1) {
		for (int i = 0; i < pos_deque.size(); i ++) {
			printf("pos_deque[%d] = %d\n", i, pos_deque[i]);
		}
	}

	// when reach the last variant, we need to calcualte the leftover geno_deque
	while (geno_deque.size() > 1) { // calcualte LD until the first in scope variant
		for (int j = 1; j < geno_deque.size(); j ++) {
			fprintf(fout, "%s\t%d\t", lastChrom.c_str(), pos_deque[0]);
			fprintf(fout, "%s\t%d\t", chrom.c_str(), pos_deque[j]);

			int n[4] = {0};
			double d[3] = {0.0}; // D prime
			if (DEBUG_Youna == 1) {
				printf("pos_deque[%d] = %d, pos_deque[%d] = %d\n", 0, pos_deque[0], j, pos_deque[j]);
			}
			calculateLD(geno_deque[0], geno_deque[j], n, d);
			fprintf(fout, "%d\t%d\t%d\t%d\t", n[0], n[1], n[2], n[3]);
			fprintf(fout, "%.6lf\t%.6lf\t%.6lf\n", d[0], d[1], d[2]);
		}
		geno_deque.pop_front();
		pos_deque.pop_front();
		if (geno_deque.size() > 0) {
			lastPos = pos_deque[0];
		}
	}

	return 1;
};
