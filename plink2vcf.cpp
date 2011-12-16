#include "Argument.h"
#include "PlinkInputFile.h"
#include "MathMatrix.h"

int main(int argc, char *argv[])
{
    time_t currentTime = time(0);
    fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

    ////////////////////////////////////////////////
    BEGIN_PARAMETER_LIST(pl)
        ADD_PARAMETER_GROUP(pl, "Input/Output")
        ADD_STRING_PARAMETER(pl, inPlink, "--inPlink", "input Plink File")
        ADD_STRING_PARAMETER(pl, outVcf, "--outVcf", "output prefix")
        // ADD_STRING_PARAMETER(pl, outPlink, "--make-bed", "output prefix")
        // ADD_PARAMETER_GROUP(pl, "People Filter")
        // ADD_STRING_PARAMETER(pl, peopleIncludeID, "--peopleIncludeID", "give IDs of people that will be included in study")
        // ADD_STRING_PARAMETER(pl, peopleIncludeFile, "--peopleIncludeFile", "from given file, set IDs of people that will be included in study")
        // ADD_STRING_PARAMETER(pl, peopleExcludeID, "--peopleExcludeID", "give IDs of people that will be included in study")
        // ADD_STRING_PARAMETER(pl, peopleExcludeFile, "--peopleExcludeFile", "from given file, set IDs of people that will be included in study")
        // ADD_PARAMETER_GROUP(pl, "Site Filter")
        // ADD_STRING_PARAMETER(pl, rangeList, "--rangeList", "Specify some ranges to use, please use chr:begin-end format.")
        // ADD_STRING_PARAMETER(pl, rangeFile, "--rangeFile", "Specify the file containing ranges, please use chr:begin-end format.")
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

    REQUIRE_STRING_PARAMETER(FLAG_inPlink, "Please provide input file using: --inVcf");

    PlinkInputFile* pin = new PlinkInputFile(FLAG_inPlink.c_str());
    FILE* fout = fopen(FLAG_outVcf.c_str(), "wt");
    if (!fout) {
        fprintf(stderr, "Cannot open output file %s!\n", FLAG_outVcf.c_str());
        abort();
    }

    int numPeople = pin->getNumIndv();
    int numMarker = pin->getNumMarker();

    fprintf(fout, "##fileformat=VCFv4.0\n");
    fprintf(fout, "##filedate=\n");
    fprintf(fout, "##source=plink2vcf\n");
    fprintf(fout, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    // writer header
    for (int p = 0; p < numPeople; p++) {
        fprintf(fout, "\t%s", pin->indv[p].c_str());
    };
    fprintf(fout, "\n");

    // write content
    Matrix mat;
    pin->readIntoMatrix(&mat);
    for (int m = 0; m < numMarker; m++){
        fprintf(fout, "%s\t", pin->chrom[m].c_str());// CHROM
        fprintf(fout, "%d\t", pin->pos[m]);          // POS
        fprintf(fout, "%s\t", pin->snp[m].c_str());  // ID
        fprintf(fout, "%c\t", pin->ref[m]);          // REF
        fprintf(fout, "%c\t", pin->alt[m]);          // ALT
        fprintf(fout, ".\t", pin->chrom[m].c_str()); // QUAL
        fprintf(fout, ".\t", pin->chrom[m].c_str()); // FILTER
        fprintf(fout, ".\t", pin->chrom[m].c_str()); // INFO
        fprintf(fout, "GT", pin->chrom[m].c_str()); // FORMAT
        for (int p = 0; p < numPeople; p++){
            int geno = mat[p][m];
            switch (geno){
            case 0:
                fputs("\t0/0", fout);
                break;
            case 1:
                fputs("\t0/1", fout);
                break;
            case 2:
                fputs("\t1/1", fout);
                break;
            default:
                fputs("\t./.", fout);
                break;
            }
        }
        fprintf(fout, "\n");
    }
    
    return 0;
}
