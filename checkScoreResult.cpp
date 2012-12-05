#include "base/Argument.h"
#include "base/GenomeSequence.h"
#include "base/IO.h"
#include "base/Logger.h"
#include "base/Utils.h"
#include "GitVersion.h"

Logger* logger = NULL;

bool reverseComplementBase(std::string* b) {
  if (b->size() != 1)
    return false;
  switch ( (*b)[0] ) {
    case 'A':
      (*b) = 'T';
      return true;
    case 'T':
      (*b) = 'A';
      return true;
    case 'G':
      (*b) = 'C';
      return true;
    case 'C':
      (*b) = 'G';
      return true;
    default:
      return false;
  }
}

/**
 * @return a string representing current time, without '\n' at the end
 */
std::string currentTime() {
  time_t t = time(NULL);
  std::string s (ctime(&t));
  s = s.substr(0, s.size() - 1);
  return s;
};

int main(int argc, char *argv[])
{
  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Input/Output")
      ADD_STRING_PARAMETER(pl, inScore, "--inScore", "input VCF File")
      ADD_STRING_PARAMETER(pl, outPrefix, "--out", "output prefix")
      ADD_STRING_PARAMETER(pl, reference, "--reference", "specify indexed reference file (.fa and .fai)")
      ADD_PARAMETER_GROUP(pl, "Auxilliary Functions")
      ADD_BOOL_PARAMETER(pl, help, "--help", "Print detailed help message")
      END_PARAMETER_LIST(pl)
      ;

  pl.Read(argc, argv);

  if (FLAG_help) {
    pl.Help();
    return 0;
  }

  pl.Status();
  if (FLAG_REMAIN_ARG.size() > 0){
    fprintf(stderr, "Unparsed arguments: ");
    for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++){
      fprintf(stderr, " %s", FLAG_REMAIN_ARG[i].c_str());
    }
    abort();
  }

  REQUIRE_STRING_PARAMETER(FLAG_inScore, "Please provide input file using: --inScore");
  REQUIRE_STRING_PARAMETER(FLAG_reference, "Please provide reference genome file using: --reference");

  if (!FLAG_outPrefix.size())
    FLAG_outPrefix = "flippedScore";

  Logger _logger( (FLAG_outPrefix + ".log").c_str());
  logger = &_logger;
  logger->infoToFile("Program Version");
  logger->infoToFile(gitVersion);
  logger->infoToFile("Parameters BEGIN");
  pl.WriteToFile(logger->getHandle());
  logger->infoToFile("Parameters END");

  time_t startTime = time(0);
  logger->info("Analysis started at: %s", currentTime().c_str());

  /* Sample score file:
   * CHROM   POS     REF     ALT     NumSample       AF      Stat    Direction       Pvalue
     1       569424  C       T       2763    NA      NA      NA      NA
     1       762320  T       C       2763    0.998371        0.0296799       -       0.863219
   */
  int totalCorrect = 0;
  int totalSwitch = 0;
  int totalDropped = 0;
  int totalFlippedCorrect = 0;
  int totalFlippedSwith = 0;  
  GenomeSequence gs;
  if (!gs.open(FLAG_reference.c_str())) {
    logger->error("Cannot open reference genome file [ %s ]", FLAG_reference.c_str());
    return -1;
  };

  std::string outFileName = FLAG_outPrefix + ".checked.SingleScore.assoc";
  FileWriter fout( outFileName.c_str() );
  
  std::string line;
  std::vector<std::string> fd;
  int lineNo = 0;
  LineReader lr(FLAG_inScore);
  while (lr.readLine(&line)) {
    lineNo ++;
    
    if (lineNo == 1) {
      if (line != "CHROM\tPOS\tREF\tALT\tNumSample\tAF\tStat\tDirection\tPvalue") {
        logger->warn("Continue anyway - the input file header is not recognized, be cautious about the results!");
      }
      fout.writeLine(line.c_str());
      continue;
    }

    if ( 9 != stringTokenize(line, "\t ", &fd)) {
      logger->error("Skip: Line [ %d ] does not have correct number of columns!", lineNo);
      totalDropped ++;      
      continue;
    };

    std::string& chrom = fd[0];
    int pos = atoi(fd[1]);
    std::string& ref = fd[2];
    std::string& alt = fd[3];
    std::string& ns = fd[4];
    std::string& af = fd[5];
    std::string& stat = fd[6];
    std::string& direction = fd[7];
    std::string& pvalue = fd[8];
    
    if (ref.size() != 1 || alt.size() != 1) {
      logger->info("Skip: Cannot check multiple bases. Ref = [ %s ], Alt = [ %s ]", ref.c_str(), alt.c_str());
      totalDropped ++;      
      continue;
    }

    if (!gs.exists(chrom)) {
      chrom = chopChr(chrom);
    }

    if (!gs.exists(chrom)) {
      logger->info("Skip: Cannot find chromosome [ %s ] in reference genome file", chrom.c_str());
      totalDropped ++;
      continue;
    }

    const char base = gs[chrom][pos - 1];
    char r = ref[0];
    char a = alt[0];

    if (r == base) { // no need to do anything
      fout.writeLine(line.c_str());
      totalCorrect ++;
      continue;
    }

    if (a == base) { // need to switch ref and alt
      // modify af, direction
      double newAf;
      if (str2double(af, &newAf) == true) {
        af = toString(1.0 - newAf);
      }
      if (direction == "+") {
        direction = "-";
      } else if (direction == "-") {
        direction = "+";
      }
      fout.printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", chrom.c_str(), pos,
                  alt.c_str(), ref.c_str(),
                  ns.c_str(), af.c_str(),
                  stat.c_str(), direction.c_str(),
                  pvalue.c_str());
                  
      
      totalSwitch ++;
      continue;
    }

    // try switch strand
    std::string rcRef = ref;
    if (!reverseComplementBase(&rcRef)) {
      logger->info("Drop site because we cannot try reverse complement of reference [ %s ]", ref.c_str());
      totalDropped ++;
      continue;
    }
    if (rcRef[0] == alt[0]) {
      logger->info("Drop site because of strand: report ref/alt [ %s/%s ], true ref [ %c ]", ref.c_str(), alt.c_str(), base);
      totalDropped ++;
      continue;
    }

    reverseComplementBase(&ref);
    reverseComplementBase(&alt);
     r = ref[0];
     a = alt[0];
    
    if (r == base) { // no need to do anything
      fout.printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", chrom.c_str(), pos,
                  ref.c_str(), alt.c_str(),
                  ns.c_str(), af.c_str(),
                  stat.c_str(), direction.c_str(),
                  pvalue.c_str());
      totalFlippedCorrect ++;
      continue;
    }

    if (a == base) { // need to switch ref and alt
      // modify af, direction
      double newAf;
      if (str2double(af, &newAf) == true) {
        af = toString(1.0 - newAf);
      }
      if (direction == "+") {
        direction = "-";
      } else if (direction == "-") {
        direction = "+";
      }
      fout.printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", chrom.c_str(), pos,
                  alt.c_str(), ref.c_str(),
                  ns.c_str(), af.c_str(),
                  stat.c_str(), direction.c_str(),
                  pvalue.c_str());
      totalFlippedSwith ++;
      continue;
    }

    logger->info("Drop site because we do not know how to correct [ %s:%d ]", chrom.c_str(), pos);
    totalDropped ++;
    
  }
  logger->info("Total processed [ %d ] lines (including header)", lineNo);
  logger->info("Total correct sites: [ %d ]", totalCorrect);  
  logger->info("Total switched sites: [ %d ]", totalSwitch);
  logger->info("Total correct sites after flipping: [ %d ]", totalFlippedCorrect);
  logger->info("Total switched sites after flipping: [ %d ]", totalFlippedSwith);
  logger->info("Total dropped sites: [ %d ]", totalDropped);
  
  time_t endTime = time(0);
  logger->info("Analysis ends at: %s", currentTime().c_str());
  int elapsedSecond = (int) (endTime - startTime);
  logger->info("Analysis took %d seconds", elapsedSecond);
  return 0;
}
