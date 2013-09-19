#include "BCFReader.h"

#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bcf.h"
#include "kstring.h"
#include "kseq.h"

  KSTREAM_INIT(gzFile, gzread, 4096)

  typedef struct {
    gzFile fp;
    FILE *fpout;
    kstream_t *ks;
    void *refhash;
    kstring_t line;
    int max_ref;
  } vcf_t;

/**
 * Adopt some functions from vcf.c
 */
extern "C" {
  void *bcf_str2id_init();
  extern void bcf_fmt_core(const bcf_hdr_t *h, bcf1_t *b, kstring_t *s);
}
static void my_write_header(bcf_hdr_t *h);

bcf_t * my_vcf_open(const char *fn, const char *mode)
{
  bcf_t *bp;
  vcf_t *v;
  if (strchr(mode, 'b')) return bcf_open(fn, mode);
  bp = (bcf_t*)calloc(1, sizeof(bcf_t));
  v = (vcf_t*)calloc(1, sizeof(vcf_t));
  bp->is_vcf = 1;
  bp->v = v;
  v->refhash = bcf_str2id_init();
  if (strchr(mode, 'r')) {
    v->fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    v->ks = ks_init(v->fp);
  } else if (strchr(mode, 'w')) {
    // disable open stdout/external file
    // v->fpout = strcmp(fn, "-")? fopen(fn, "w") : stdout;
    v->fpout = NULL;
  }
  return bp;
}

#define my_vcf_read vcf_read

// vcf_hdr_read
// vcf_hdr_write(bout, hout, &header);

int my_vcf_hdr_write(bcf_t *bp, const bcf_hdr_t *h, std::string* hdr) {
  // vcf_t *v = (vcf_t*)bp->v;
  int i, has_ver = 0;
  if (!bp->is_vcf) {
    fprintf(stderr, "Something is wrong when reading BCF header at %s:%d\n", __FILE__, __LINE__);
    return bcf_hdr_write(bp, h);
  }
  std::string& s = *hdr;
  if (h->l_txt > 0) {
    if (strstr(h->txt, "##fileformat=")) has_ver = 1;
    if (has_ver == 0) {
      //fprintf(v->fpout, "##fileformat=VCFv4.1\n");
      s = "##fileformat=VCFv4.1\n";
    }
    // fwrite(h->txt, 1, h->l_txt - 1, v->fpout);
    s += h->txt;
  }
  if (h->l_txt == 0) {
    //fprintf(v->fpout, "##fileformat=VCFv4.1\n");
    s = "##fileformat=VCFv4.1\n";
  }
  // fprintf(v->fpout, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
  s += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for (i = 0; i < h->n_smpl; ++i) {
    //fprintf(v->fpout, "\t%s", h->sns[i]);
    s += "\t";
    s += h->sns[i];
  }
  //fputc('\n', v->fpout);
  s += "\n";
  return 0;
}

// adopted from vcf_write() in vcf.c
int my_vcf_write(bcf_t *bp, bcf_hdr_t *h, bcf1_t *b, std::string* line) {
  // vcf_t *v = (vcf_t*)bp->v;

  if (!bp->is_vcf) {
    fprintf(stderr, "Something is wrong when reading BCF at %s:%d\n", __FILE__, __LINE__);
    return bcf_write(bp, h, b);
  }

  kstring_t str;
  bcf_fmt_core(h, b, &str);
  // bcf_fmt_core(h, b, &v->line);
  // fwrite(v->line.s, 1, v->line.l, v->fpout);
  // fputc('\n', v->fpout);
  line->assign(str.s, str.l);
  return str.l + 1;
  // return v->line.l + 1;
}

int BCFReader::open(const std::string& fn) {
  inReading = false;

  // check file existance
  bp = my_vcf_open(fn.c_str(), "rb"); // vcf file handle
  if (!bp) {
    this->cannotOpen = true;
    return -1;
  }
  b = (bcf1_t*) calloc(1, sizeof(bcf1_t)); // each bcf read block
  if (b == 0) {
    return -1;
  }
  // write header
  hin = hout = vcf_hdr_read(bp);

  bout = my_vcf_open("-", "wu");

  my_write_header(hout); // always print the header, put certain fields in header.

  // write results out
  // vcf_hdr_write(bout, hout);
  my_vcf_hdr_write(bout, hout, &header);

  //open index
  this->hasIndex = this->openIndex(fn);

  // set up range iterator
  resetRangeIterator();

  cannotOpen = false;
  return 0;
}


bool BCFReader::readLine(std::string* line) {
  // openOK?
  if (cannotOpen) return false;

  // read
  if (!inReading) {
    resetRangeIterator();
    inReading = true;
  };

  // check read mode
  if (range.empty()) {
    // read line by line
    if (my_vcf_read(bp, hin, b) > 0) {
      my_vcf_write(bout, hout, b, line);
      return true;
    }
    return false;
  }

  // read by region
  // check index
  assert(!range.empty());
  if (!hasIndex) return false;

  if (off != 0) {
    // read a record
    while (my_vcf_read(bp, hin, b) > 0) {
      if (tid >= 0) {
        int l = strlen(b->ref);
        l = b->pos + (l > 0? l : 1);
        if (b->tid != tid || b->pos >= end) break; // current record has passed prespecified region
        if (!(l > begin && end > b->pos)) continue; // not sure when this will happen

        my_vcf_write(bout, hout, b, line);
        return true;
      }
    }
  }

  // seek to region
  for (; this->rangeIterator != this->rangeEnd; ++ rangeIterator) {
    char rangeBuffer[128];
    snprintf(rangeBuffer, 128, "%s:%u-%u", this->rangeIterator.getChrom().c_str(),
             this->rangeIterator.getBegin(), this->rangeIterator.getEnd());
    rangeBuffer[127] = '\0';
    // int tid, beg, end, len;
    if (!str2id) {
      str2id = bcf_build_refhash(hout);
    }
    if (bcf_parse_region(str2id, rangeBuffer, &tid, &begin, &end) >= 0) {
      off = bcf_idx_query(idx, tid, begin);
      if (off == 0) {
        // fprintf(stderr, "[%s] no records in the query region.\n", __func__);
        // return 1; // FIXME: a lot of memory leaks...
        continue;
      }
      bgzf_seek(bp->fp, off, SEEK_SET);

      while (my_vcf_read(bp, hin, b) > 0) {
        if (tid >= 0) {
          int l = strlen(b->ref);
          l = b->pos + (l > 0? l : 1);
          if (b->tid != tid || b->pos >= end) break; // current record has passed prespecified region
          if (!(l > begin && end > b->pos)) continue; // not sure when this will happen

          ++rangeIterator;
          my_vcf_write(bout, hout, b, line);
          return true;
        }
      }
    }
  }
  return false;
};

static void my_write_header(bcf_hdr_t *h)
{
  kstring_t str;
  str.l = h->l_txt? h->l_txt - 1 : 0;
  str.m = str.l + 1;
  str.s = h->txt;
  if (!strstr(str.s, "##INFO=<ID=DP,"))
    kputs("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=DP4,"))
    kputs("##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=MQ,"))
    kputs("##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root-mean-square mapping quality of covering reads\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=FQ,"))
    kputs("##INFO=<ID=FQ,Number=1,Type=Float,Description=\"Phred probability of all samples being the same\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=AF1,"))
    kputs("##INFO=<ID=AF1,Number=1,Type=Float,Description=\"Max-likelihood estimate of the first ALT allele frequency (assuming HWE)\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=AC1,"))
    kputs("##INFO=<ID=AC1,Number=1,Type=Float,Description=\"Max-likelihood estimate of the first ALT allele count (no HWE assumption)\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=AN,"))
    kputs("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=IS,"))
    kputs("##INFO=<ID=IS,Number=2,Type=Float,Description=\"Maximum number of reads supporting an indel and fraction of indel reads\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=AC,"))
    kputs("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes for each ALT allele, in the same order as listed\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=G3,"))
    kputs("##INFO=<ID=G3,Number=3,Type=Float,Description=\"ML estimate of genotype frequencies\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=HWE,"))
    kputs("##INFO=<ID=HWE,Number=1,Type=Float,Description=\"Chi^2 based HWE test P-value based on G3\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=CLR,"))
    kputs("##INFO=<ID=CLR,Number=1,Type=Integer,Description=\"Log ratio of genotype likelihoods with and without the constraint\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=UGT,"))
    kputs("##INFO=<ID=UGT,Number=1,Type=String,Description=\"The most probable unconstrained genotype configuration in the trio\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=CGT,"))
    kputs("##INFO=<ID=CGT,Number=1,Type=String,Description=\"The most probable constrained genotype configuration in the trio\">\n", &str);
  //	if (!strstr(str.s, "##INFO=<ID=CI95,"))
  //		kputs("##INFO=<ID=CI95,Number=2,Type=Float,Description=\"Equal-tail Bayesian credible interval of the site allele frequency at the 95% level\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=PV4,"))
    kputs("##INFO=<ID=PV4,Number=4,Type=Float,Description=\"P-values for strand bias, baseQ bias, mapQ bias and tail distance bias\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=INDEL,"))
    kputs("##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=PC2,"))
    kputs("##INFO=<ID=PC2,Number=2,Type=Integer,Description=\"Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than in group2.\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=PCHI2,"))
    kputs("##INFO=<ID=PCHI2,Number=1,Type=Float,Description=\"Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples.\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=QCHI2,"))
    kputs("##INFO=<ID=QCHI2,Number=1,Type=Integer,Description=\"Phred scaled PCHI2.\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=RP,"))
    kputs("##INFO=<ID=PR,Number=1,Type=Integer,Description=\"# permutations yielding a smaller PCHI2.\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=QBD,"))
    kputs("##INFO=<ID=QBD,Number=1,Type=Float,Description=\"Quality by Depth: QUAL/#reads\">\n", &str);
  //if (!strstr(str.s, "##INFO=<ID=RPS,"))
  //    kputs("##INFO=<ID=RPS,Number=3,Type=Float,Description=\"Read Position Stats: depth, average, stddev\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=RPB,"))
    kputs("##INFO=<ID=RPB,Number=1,Type=Float,Description=\"Read Position Bias\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=MDV,"))
    kputs("##INFO=<ID=MDV,Number=1,Type=Integer,Description=\"Maximum number of high-quality nonRef reads in samples\">\n", &str);
  if (!strstr(str.s, "##INFO=<ID=VDB,"))
    kputs("##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Variant Distance Bias (v2) for filtering splice-site artefacts in RNA-seq data. Note: this version may be broken.\">\n", &str);
  if (!strstr(str.s, "##FORMAT=<ID=GT,"))
    kputs("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n", &str);
  if (!strstr(str.s, "##FORMAT=<ID=GQ,"))
    kputs("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n", &str);
  if (!strstr(str.s, "##FORMAT=<ID=GL,"))
    kputs("##FORMAT=<ID=GL,Number=3,Type=Float,Description=\"Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)\">\n", &str);
  if (!strstr(str.s, "##FORMAT=<ID=DP,"))
    kputs("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"# high-quality bases\">\n", &str);
  if (!strstr(str.s, "##FORMAT=<ID=DV,"))
    kputs("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality non-reference bases\">\n", &str);
  if (!strstr(str.s, "##FORMAT=<ID=SP,"))
    kputs("##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">\n", &str);
  if (!strstr(str.s, "##FORMAT=<ID=PL,"))
    kputs("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n", &str);
  h->l_txt = str.l + 1;
  h->txt = str.s;
}
