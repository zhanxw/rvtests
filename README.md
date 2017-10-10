<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [Introduction](#introduction)
- [Citation](#citation)
- [Download](#download)
- [Quick tutorial](#quick-tutorial)
    - [Single variant tests](#single-variant-tests)
    - [Groupwise tests](#groupwise-tests)
    - [Related individual tests](#related-individual-tests)
    - [Meta-analysis tests](#meta-analysis-tests)
        - [Dominant models and recessive models](#dominant-models-and-recessive-models)
- [Input files](#input-files)
    - [Genotype files (VCF, BCF, BGEN, KGG)](#genotype-files-vcf-bcf-bgen-kgg)
    - [Phenotype file](#phenotype-file)
    - [Covariate file](#covariate-file)
    - [Trait transformation](#trait-transformation)
- [Models](#models)
    - [Single variant tests](#single-variant-tests-1)
    - [Burden tests](#burden-tests)
    - [Variable threshold models](#variable-threshold-models)
    - [Kernel models](#kernel-models)
    - [Meta-analysis models](#meta-analysis-models)
    - [Utility models](#utility-models)
- [Association test options](#association-test-options)
    - [Sample inclusion/exclusion](#sample-inclusionexclusion)
    - [Variant site filters](#variant-site-filters)
        - [Annotation](#annotation)
    - [Genotype filters](#genotype-filters)
    - [Handle missing genotypes and phenotypes](#handle-missing-genotypes-and-phenotypes)
    - [Specify groups (e.g burden unit)](#specify-groups-eg-burden-unit)
- [Sex chromosome analysis](#sex-chromosome-analysis)
- [Kinship generation](#kinship-generation)
- [Resources](#resources)
    - [UCSC RefFlat Genes](#ucsc-refflat-genes)
    - [Gencode Genes](#gencode-genes)
- [Frequently Asked Questions (FAQ)](#frequently-asked-questions-faq)
- [Feedback/Contact](#feedbackcontact)

<!-- markdown-toc end -->

[![Build Status](https://travis-ci.org/zhanxw/rvtests.png?branch=master)](https://travis-ci.org/zhanxw/rvtests)

(Updated: October 2017)

# Introduction

Rvtests, which stands for Rare Variant tests, is a flexible software package for genetic association analysis for sequence datasets. Since its inception, rvtests was developed as a comprehensive tool to support genetic association analysis and meta-analysis. It can analyze both unrelated individual and related (family-based) individuals for both quantitative and binary outcomes. It includes a variety of association tests (e.g. single variant score test, burden test, variable threshold test, SKAT test, fast linear mixed model score test). It takes [VCF][vcf]/BGEN/PLINK format as genotype input file and takes PLINK format phenotype file and covariate file. 

With new implementation of the BOLT-LMM/MINQUE algorithm as well as a series of software engineering optimizations, our software package is capable of analyzing datasets of up to 1,000,000 individuals in linear mixed models on a computer workstation, which makes our tool one of the very few options for analyzing large biobank scale datasets, such as UK Biobank. RVTESTS supports both single variant and gene-level tests. It also allows for highly effcient generation of covariance matrices between score statistics in RAREMETAL format, which can be used to support the next wave of meta-analysis that incorporates large biobank datasets.

A (much) larger sample size can be handled using linear regression or logistic regression models. 

[vcf]: http://www.1000genomes.com/

# Citation

Xiaowei Zhan, Youna Hu, Bingshan Li, Goncalo R. Abecasis, and Dajiang J. Liu

**RVTESTS: An Efficient and Comprehensive Tool for Rare Variant Association Analysis Using Sequence Data**

Bioinformatics 2016 32: 1423-1426. [doi:10.1093/bioinformatics/btw079](http://bioinformatics.oxfordjournals.org/content/32/9/1423.short)  ([PDF](http://bioinformatics.oxfordjournals.org/content/32/9/1423.full.pdf+html))

# Download

Source codes can be downloaded from [github](https://github.com/zhanxw/rvtests/archive/master.zip) or [github page](https://github.com/zhanxw/rvtests). In Linux, you can use 
```
git clone https://github.com/zhanxw/rvtests
``` 
to retrieve the latest distribution for rvtests. To install, go to the rvtests folder and type `make`. When compilation succeed, the executable is under the `executable` folder. Simply type `executable/rvtests` can get you started.

Alternatively, binary executable files (for Linux 64-bit platform) can be downloaded from [here](https://github.com/zhanxw/rvtests/releases).


# Quick tutorial

Here is a quick example of how to use *rvtests* software in typical use cases.

## Single variant tests

    rvtest --inVcf input.vcf --pheno phenotype.ped --out output --single wald,score

This specifies single variant Wald and score test for association
tests for every variant in the `input.vcf` file. The 6th column of the phenotype file, `phenotype.ped`, which is in PLINK format, is used. Rvtests will automatically check whether the phenotype is binary trait or quantitative trait.
For binary trait, the recommended way of coding is to code controls as 1, cases as 2, missing phenotypes as -9 or 0.

For other types of association tests, you can refer to [Models](#models).

## Groupwise tests
Groupwise tests includes three major kinds of tests.

* Burden tests: group variants, which are usually less than 1% or 5% rare variants, for association tests. The category includes: CMC test, Zeggini test, Madsen-Browning test, CMAT test, and rare-cover test.
* Variable threshold tests: group variants under different frequency thresholds.
* Kernel methods: suitable to tests rare variants having different directions of effects. These includes SKAT test and KBAC test. 

All above tests requires to group variants into a unit. The simplest case is to use gene as grouping unit. For different grouping method, see [Grouping](#Grouping). 

To perform rare variant tests by gene, you need to use `--geneFile` to specify the gene range in a refFlat format. We provided different gene definitions in the [Resources](#Resources) section. You can use `--gene` to specify which gene(s) to test. For example, specify `--gene CFH,ARMS2` will perform association tests on CFH and ARMS2 genes. If there is no providing `--gene` option, all genes will be tests.

The following command line demonstrate how to use CMC method, variable threshold method(proposed by Price) and kernel based method (SKAT by Shawn Lee and KBAC by
Dajiang Liu) to test every gene listed in *refFlat\_hg19.txt.gz*.

    rvtest --inVcf input.vcf --pheno phenotype.ped --out output --geneFile refFlat_hg19.txt.gz --burden cmc --vt price --kernel skat,kbac


## Related individual tests

To test related individuals, you will need to first create a kinship matrix:

    vcf2kinship --inVcf input.vcf --bn --out output

The option `--bn` means calculating empirical kinship using Balding-Nicols method. You can specify `--ibs` to obtain IBS kinship or use `--pedigree input.ped` to calculate kinship from known pedigree information.

Then you can use linear mixed model based association tests such as Fast-LMM score test, Fast-LMM LRT test and Grammar-gamma tests. An exemplar command is shown: 

    rvtest --inVcf input.vcf --pheno phenotype.ped --out output --kinship output.kinship --single famScore,famLRT,famGrammarGamma

## Meta-analysis tests

The meta-analysis models outputs association test results and genotype covariance matrix. 
These calculated summary statistics can be used in rare variant association analysis ([details](#meta-analysis-models)).
We provide single variant score test and generate a genotype covariance matrix. 
You can use this command:
   
    rvtest --inVcf input.vcf --pheno phenotype.ped --meta score,cov --out output

In a more realistic scenario, you may want to adjust for covariates and want to inverse normalized residuals obtained in null model ([link](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.2852.html) to our methodology paper), then this command will work:
   
    rvtest --inVcf input.vcf --pheno phenotype.ped --covar example.covar --covar-name age,bmi --inverseNormal --useResidualAsPhenotype  --meta score,cov --out output


Here the `--covar` specify a covariate file, and `--covar-name` specify which covariates can used in the analysis. Covariate file format can be found [here](#Covariate file). `--inverseNormal --useResidualAsPhenotype` specifies trait transformation method. That means first fit a regression model of the phenotype on covariates (intercept automatically added), then the residuals are inverse normalized. Trait transformation details can be found [here](#Trait transformation).

We support both unrelated individuals and related individuals (e.g. family data). You need to append `--kinship input.kinship` to the command line:

    rvtest --inVcf input.vcf --pheno phenotype.ped --meta score,cov --out output --kinship input.kinship

The file `input.kinship` is calculated by `vcf2kinship` program, and usage to this program is described in [Related individual tests](#related-individual-tests).

**NOTE:** by default, the covariance matrix are calculated in a sliding-window of 1 million base pairs. You can change this setting via  the option `windowSize`.
For example, `--meta cov[windowSize=500000]` specify a 500k-bp sliding window.


### Dominant models and recessive models

Dominant and recessive disease models are supported by appending "dominant" and/or "recessive" after "--meta" option. For example, use "--meta dominant,recessive" will 
generate two sets of files. For dominant model, they are "prefix.MetaDominant.assoc" and "prefix.MetaDominantCov.assoc.gz"; for recessive model,
they are "prefix.MetaRecessive.assoc" and "prefix.MetaRecessiveCov.assoc.gz". Internally, in dominant models, genotypes 0/1/2 are coded as 0/1/1; in recessive models, genotypes 0/1/2 are 
coded as 0/0/1. Missing genotypes will be imputed to the mean.

# Input files

## Genotype files (VCF, BCF, BGEN, KGG)

Rvtests supports VCF (Variant Call Format) files. Files in both plain text format or gzipped format are supported. To use group-based rare variant tests, indexed the VCF files using [tabix](http://samtools.sourceforge.net/tabix.shtml) are required. 

Here are the commands to convert plain text format to bgzipped VCF format:

    (grep ^"#" $your_old_vcf; grep -v ^"#" $your_old_vcf | sed 's:^chr::ig' | sort -k1,1n -k2,2n) | bgzip -c > $your_vcf_file 
    tabix -f -p vcf $your_vcf_file

The above commands will (1) remove the `chr` prefix from chromosome names; (2) sort VCF files by chromosome first, then by chromosomal positions; (3) compress using bgzip; (4) create tabix index.

Rvtests support genotype dosages. Use `--dosage DosageTag` to specify the dosage tag. For example, if VCF format field is "GT:EC" and individual genotype fields is "0/0:0.02", you can use `--dosage EC`, and rvtests will use the dosage 0.02 in the regression models.

Rvtests suppport [BGEN input format](http://www.well.ox.ac.uk/~gav/bgen_format/) v1.0 throught v1.3. Instead of using `--inVcf`, use `--inBgen` to specify a BGEN file and `--inBgenSample` to specify the accompany SAMPLE file.

Rvtests support [KGGSeq input format](http://grass.cgs.hku.hk/limx/kggseq/doc10/UserManual.html#kedformat). This format is an extension to binary PLINK formats. Use `--inKgg` to replace `--inVcf`.


## Phenotype file

You can use `--mpheno $phenotypeColumnNumber` or `--pheno-name` to specify a given phenotype.

An example phenotype file, (`example.pheno`), has the following format: 

    fid iid fatid matid sex y1 y2 y3 y4
    P1 P1 0 0 0 1.7642934435605 -0.733862638327895 -0.980843608339726 2
    P2 P2 0 0 0 0.457111744989746 0.623297281416372 -2.24266162284447 1
    P3 P3 0 0 0 0.566689682543218 1.44136462889459 -1.6490100777089 1
    P4 P4 0 0 0 0.350528353203767 -1.79533911725537 -1.11916876241804 1
    P5 P5 0 0 1 2.72675074738545 -1.05487747371158 -0.33586430010589 2

Phenotype file is specified by the option `--pheno example.pheno` . The default phenotype column header is “`y1`”. If you want to use alternative columns as phenotype for association analysis (e.g the column with header y2), you may specify the phenotype by column or by name using either

* --mpheno 2 
* --pheno-name y2

**NOTE:** to use “`--pheno-name`”, the  header line must starts with “`fid iid`” as PLINK requires.

In phenotype file, missing values can be denoted by NA or any non-numeric values. Individuals with missing phenotypes will be automatically dropped from subsequent association analysis. For each missing phenotype value, a warning will be generated and recorded in the log file.

When the phenotype values are only 0, 1 and 2, rvtests will automatically treat it as binary traits. However, if you want to treat it as continuous trait, please use "`--qtl`" option.

## Covariate file

You can use `--covar` and `--covar-name` to specify covariates that will be used for single variant association analysis. This is an optional parameter. If you do not have covariate in the data, this option can be ignored. 

The covariate file, (e.g. `example.covar`) has a similar format as the phenotype file:

    fid iid fatid matid sex y1 y2 y3 y4
    P1 P1 0 0 0 1.911 -1.465 -0.817 1
    P2 P2 0 0 0 2.146 -2.451 -0.178 2
    P3 P3 0 0 0 1.086 -1.194 -0.899 1
    P4 P4 0 0 0 0.704 -1.052 -0.237 1
    P5 P5 0 0 1 2.512 -3.085 -2.579 1

The covariate file is specified by the `--covar` option (e.g. `--covar example.covar`). To specify covariates that will be used in the association analysis, the option `--covar-name` can be used. For example, when age, bmi and 3 PCs are used for association analysis, the following option can be specified for the rvtests program, i.e. 
`--covar example.covar --covar-name age,bmi,pc1,pc2,pc3`.

Note: Missing data in the covariate file can be labeled by any non-numeric value (e.g. NA). They will be automatically imputed to the mean value in the data file. 


## Trait transformation

In this meta-analysis, we use inverse normal transformed residuals in the association analysis, which is achieved by using a combination of `--inverseNormal`  and `--useResidualAsPhenotype`. Specifically, we first fit the null model by regressing phenotype on covariates. The residuals are then inverse normal transformed (see Appendix A more detailed formula for transformation). Transformed residuals will be used to obtain score statistics. 

In meta analysis, an exemplar command for using rvtests looks like the following:

    ./rvtest --inVcf $vcf --pheno $example.pheno --covar example.covar --covar-name age,bmi --inverseNormal --useResidualAsPhenotype  --meta score,cov --out $output_prefix  

# Models

Rvtests support various association models.

## Single variant tests

Single variant | Model(#)    |Traits(##) | Covariates | Related / unrelated | Description
:--------------|:---------:|:------:|:----------:|:-------------------:|:-----------
Score test     |  score    |B, Q  |     Y      |         U           | Only null model is used to performed the test
Wald  test     |  wald     |B, Q  |     Y      |         U           | Only fit alternative model, and effect size will be estimated
Exact test     |  exact    |B     |     N      |         U           | Fisher's test (allelic test)
Dominant Exact test     |  dominantExact    |B     |     N      |         U           | Fisher's test (dominant codings)
Fam LRT        |  famLRT   |Q     |     Y      |         R, U        | Fast-LMM model
Fam Score      |  famScore |Q     |     Y      |         R, U        | Fast-LMM model style likelihood ratio test
Grammar-gamma  |famGrammarGamma| Q     |     Y      |         R, U        | Grammar-gamma method
Firth regression  |firth| B     |     Y      |         U        | Logistic regression with Firth correction by David Firth, discussed by Clement Ma.

(#) Model columns list the recognized names in rvtests. For example, use `--single score` will apply score test.

(##) In trait column, B or Q stand for binary or quantitative trait, respectively.


## Burden tests

Burden tests | Model(#)    |Traits(##) | Covariates | Related / unrelated | Description
:--------------|:---------:|:------:|:----------:|:-------------------:|:-----------
CMC             |  cmc       |B, Q  |     Y      |         U           | Collapsing and combine rare variants by Bingshan Li.
Zeggini         |  zeggini   |B, Q  |     Y      |         U           | Aggregate counts of rare variants by Morris Zeggini.
Madsen-Browning |  mb        |B     |     Y      |         U           | Up-weight rare variant using inverse frequency from controls by Madsen.
Fp              |  fp        |B     |     Y      |         U           | Up-weight rare variant using inverse frequency from controls by Danyu Lin.
Exact CMC       |  exactCMC  |B     |     N      |         U           | Collapsing and combine rare variants, then perform Fisher's exact test.
CMC Wald        |  cmcWald   |B, Q  |     Y      |         U           | Collapsing and combine rare variants, then perform Wald test.
RareCover       |  rarecover |B     |     N      |         U           | Find optimal grouping unit for rare variant tests by Thomas Hoffman.
CMAT            |  cmat      |B     |     N      |         U           | Test non-coding variants by Matt Zawistowski.
FamCMC          |  famcmc       |Q  |     Y      |         R           | Collapsing and combine rare variants extended to related samples.
FamZeggini      |  famzeggini   |Q  |     Y      |         R           | Aggregate counts of rare variants extended to related samples.


(#) Model columns list the recognized names in rvtests. For example, use `--burden cmc` will apply CMC test.

(##) In trait column, B or Q stand for binary or quantitative trait, respectively.


## Variable threshold models

Single variant | Model(#)    |Traits(##) | Covariates | Related / unrelated | Description
:--------------|:---------:|:------:|:----------:|:-------------------:|:-----------
Variable threshold model by permutation     |  price    |B, Q  |     N      |         U           | Every rare-variant frequency cutoffs are tests by Alkes Price.
Variable threshold model by analytic form   |  analytic    | Q  |     Y      |         U   | Every rare-variant frequency cutoffs are tests by Danyu Lin.
Variable threshold model by analytic form   |  famAnalytic    | Q  |     Y      |         R   | Every rare-variant frequency cutoffs are tests by Dajiang Liu.

(#) Model columns list the recognized names in rvtests. For example, use `--vt price` will variable threshold test.
NOTE: our implementation of Price's test diffs from the original method described in Price's publication. We test every minor allele frequency cutoff (instead of reference allele counts) and this is a two-sided test (instead of one-sided test).

(##) In trait column, B or Q stand for binary or quantitative trait, respectively.


## Kernel models

Kernel | Model(#)    |Traits(##) | Covariates | Related / unrelated | Description
:--------------|:---------:|:------:|:----------:|:-------------------:|:-----------
SKAT     |  skat    |B, Q  |     Y      |         U           | Sequencing kernel association test by Shawn Lee.
SKATO     |  skato    |B, Q  |     Y      |         U           | Optimal sequencing kernel association test (SKAT-O) by Shawn Lee. (###)
KBAC     |  kbac     |B  |     N      |         U           | Kernel-based adaptive clustering model by Dajiang Liu.
FamSKAT     |  famSkat    |Q  |     Y      |         R           | Sequencing kernel association test extended to related individuals by Han Chen.

(#) Model columns list the recognized names in rvtests. For example, use `--kernel skat` will apply SKAT test.
To further customize SKAT test, you can use *--kernel skat[nPerm=100:alpha=0.001:beta1=1:beta2=20]* to specify permutation counts, type-1 error, 
beta distribution parameters for up-weighting rare variants. Rvtests will output a message showing: 

    [INFO]  SKAT test significance will be evaluated using 10000 permutations at alpha = 0.001 (beta1 = 1.00, beta2 = 20.00)

Tip: When specifying beta1=1 and beta2 = 1, the associate test is equivalent to an unweighted SKAT test.

(##) In trait column, B or Q stand for binary or quantitative trait, respectively.

(###) SKAT-O implementation may have slightly different results compared with outputs from SKAT R package. That's probably due to numerical stability.

## Meta-analysis models

Type | Model(#)    |Traits(##) | Covariates | Related / unrelated | Description
:--------------|:---------:|:------:|:----------:|:-------------------:|:-----------
Score test           |  score      | B,Q  |     Y      |         R, U           | standard score tests
Dominant model       |  dominant   | B,Q  |     Y      |         R, U           | score tests and covariance matrix under dominant disease model
Recessive model      |  recessive  | B,Q  |     Y      |         R, U           | score tests and covariance matrix under recessive disease model
Covariance           |  cov        | B,Q  |     Y      |         R, U           | covariance matrix
BOLT-LMM score test           |  bolt      | Q  |     Y      |         R           | BOLT-LMM based score tests (###)
BOLT-LMM covariance           |  boltCov      | Q  |     Y      |         R           | BOLT-LMM based score tests (###)

(#) Model columns list the recognized names in rvtests. For example, use `--meta score,cov` will generate score statistics and covariance matrix for meta-analysis.

(##) In trait column, B or Q stand for binary or quantitative trait, respectively.

(###) This is an experimental feature. 
This method requires LD-pruned genotype data in the binary PLINK format specified by `--boltPlink`. 
A minimal example to run BOLT-LMM based score tests is: `rvtest --inVcf $inputVCFfile --boltPlink $binaryPlinkPrefix --pheno $phenotype --meta bolt`. 
To prune your genotype data, an example command line is `plink --vcf $inputVCFFile --maf 0.05 --indep-pairwise 50 5 0.5 --make-bed $binaryPlinkPrefix`.
Please note RVTESTS additionally prohibit duplicated rs IDs in the PLINK BIM file. 
To remove duplications, you can list all duplicated `cut -f2 $binaryPlinkBIMfile|sort |uniq -d > duplicate.rsid`, and then 
use `plink --vcf $inputVCFFile --maf 0.05 --indep-pairwise 50 5 0.5 --exclude duplicate.rsid --make-bed $binaryPlinkPrefix`.
An example script to run BOLT-LMM model can be found under `example/experimental` directory.

The above models are suitable to generate summary statistics which can be later meta-analyzed (see [Dajiang Liu (2014) Nature Genetics](http://www.nature.com/ng/journal/v46/n2/abs/ng.2852.html)).
Rvtests implemented the above methods and the results can be further analyzed by RareMetals ([link](http://genome.sph.umich.edu/wiki/RareMETALS)) for quantitative trait and RareMetals2 ([link](http://genome.sph.umich.edu/wiki/RareMETALS2)).
It also worth to mention that our group offers another toolset for meta analysis ([link](http://genome.sph.umich.edu/wiki/Rare-Metal)).

**Explanation of outputs**

The meta-analysis results come as two files: summary score statistics files (prefix.MetaScore.assoc.gz) and covariance files (prefix.MetaCov.assoc.gz).

In summary score statistics files, you will obtain these columns:

- N_INFORMATIVE: Number of samples that are analyzed for association. 
- AF: allele frequency. For related individuals, we use BLUE estimator. For case-control study, we list overall frequency (adjusted by relatedness if possible), case frequency and control frequency separated by a colon.
- INFORMATIVE_ALT_AC: The number of alternative alleles in the analyzed samples.
- CALL_RATE: The fraction of non-missing alleles. For case-control study, we calculate call rate for all samples, case samples and control samples separated by a colon.
- HWE_PVALUE: Hardy-Weinberg equilibrium. For related individuals, this statistic can be inflated. For case-control study, we calculate HWE pvalues for all samples, case samples and controls samples separated by a colon.
- N_REF/N_HET/N_ALT: Number of samples carrying homozygous reference/heterozygous/homozygous alternative alleles. For case-control study, we calculate these three statistics for all samples, case samples and controls samples separated by a colon.
- U_STAT, SQRT_V_STAT: U and V statistics are score statistics and their covariance matrix. Details can be found in [Dajiang Liu (2014) Nature Genetics](http://www.nature.com/ng/journal/v46/n2/abs/ng.2852.html).
- ALT_EFFSIZE: for continuous outcome, this is the estimated effect size; for binary outcome, this is the estimated log odds-ratio. We apply a new correction method when binary trait associations for
related individuals are analyzed in standard linear mixed models. The log odds ratio is approximately correct for related individual analysis as well.

The covariance files stores covariances in sliding windows. You will obtain these columns:

- CHROM: the chromosome name
- START_POS/END_POS: within one sliding window, the first/last variant chromosomal position.
- NUM_MARKER: number of markers in one sliding window.
- MARKER_POS: all chromosomal positions in the sliding window.
- COV: covariances between markers. When the sliding window has markers G.1, G.2, ..., G.N, this column will store "covariances" of (G.1, G.1), (G.1, G.2), ..., (G.1, G.N);
for binary outcomes, this column uses colons as separators and additionally stores covariances between N genotypes and K covariates (G.1, C.1), (G.1, C.2), ..., (G.1, C.K), as well as
covariances between K covariates (C.1, C.1), (C.1, C.2), ... (C.1, C.K), (C.2, C.2), ..., (C.2, C.K), ..., (C.K, C.K).


## Utility models


Rvtests has an convenient option `--outputRaw`. When specifying this, rvtests can output genotypes, phenotype, covariates (if any) and collapsed genotype to tabular files. These files can be imported into other software (e.g. R) for further analyses.


# Association test options

## Sample inclusion/exclusion

Rvtests can flexibly specify which sample(s) to include or exclude:

           --peopleIncludeID : give IDs of people that will be included in study
         --peopleIncludeFile : from given file, set IDs of people that will be included in study
           --peopleExcludeID : give IDs of people that will be included in study
         --peopleExcludeFile : from given file, set IDs of people that will be included in study

`--peopleIncludeID` and `--peopleExcludeID` are used to include/exclude samples from command line. 
For example, specify `--peopleIncludeID A,B,C` will include A, B and C sample from the VCF files if they exists.
`--peopleIncludeID` and `--peopleExcludeID` followed by a file name will include or exclude the IDs in the file.
So to include sample A, B and C, you can provide a file, `people.txt`, which looks like:

    A
    B
    C

Then use `--peopleIncludeFile people.txt` to include them in the analysis.


## Variant site filters

It is common that different frequency cutoffs are applied in rare-variant analysis.
Therefore, rvtests specify frequency cutoffs.

Frequency Cutoff

                 --freqUpper : Specify upper minor allele frequency bound to be included in analysis
                 --freqLower : Specify lower minor allele frequency bound to be included in analysis

If you specify `--freqLower 0.01 --freqUpper 0.05`, only the variants with minor allele frequency between 0.01 and 0.05 (boundary inclusive) will be analyzed.

Similar to sample inclusion/exclusion options, you can specify a range of variants to be included by 
specifying `--rangeList` option. For example `--rangeList 1:100-200` will include the chromosome 1 position 100bp to 200bp region.
Alternatively, use a separate file, `range.txt`, and `--rangeFile range.txt` to specify association tests range.

                 --rangeList : Specify some ranges to use, please use chr:begin-end format.
                 --rangeFile : Specify the file containing ranges, please use chr:begin-end format.
                  --siteFile : Specify the file containing sites to include, please use "chr pos" format.

It is supported to filter variant site by site depth, minor allele count or annotation (annotated VCF file is needed).

              --siteDepthMin : Specify minimum depth(inclusive) to be included in analysis
              --siteDepthMax : Specify maximum depth(inclusive) to be included in analysis
                --siteMACMin : Specify minimum Minor Allele Count(inclusive) to be included in analysis
                  --annoType : Specify annotation type that is followed by ANNO= in the VCF INFO field, regular expression is allowed

*NOTE*: `--annoType` can filter variants based on regular expression.
For example, `--annoType Nonsynonymous` will only analyze non-synonymous variants where they have `ANNO=Nonsynonymous` in the INFO field.
To extract more than one annotation types, use `--annoType 'Start_Gain|Stop_Loss|Start_Loss|Essential_Splice_Site|Stop_Gain|Normal_Splice_Site|Synonymous|Nonsynonymous` will extract LOF (loss of function) mutations.
To generate an annotated VCF file, please read on the next section [Annotation](#Annotation). 


### Annotation

We define a VCF file with annotating information as an annotated VCF here.
The annotation step can be done using [ANNO](https://github.com/zhanxw/anno), a fast and accurate annotation software.
At minimal, there are two steps to annotate a VCF file:

1. Install ANNO and its resources files

    git clone https://github.com/zhanxw/anno.git;
    cd anno; make;
    cd resources; ./download.sh

2. Run the following script:

    anno -i input.vcf -o output.vcf.gz -r anno/resources/hs37d5.fa -g anno/resources/refFlat_hg19.txt.gz -p anno/priority.txt -c anno/codon.txt --indexOutput

You will then obtain an annotated VCF file `output.vcf.gz` and its tabix index `output.vcf.gz.tbi`.


## Genotype filters

Genotype with low depth or low quality can be filtered out by these options:

              --indvDepthMin : Specify minimum depth(inclusive) of a sample to be included in analysis
              --indvDepthMax : Specify maximum depth(inclusive) of a sample to be included in analysis
               --indvQualMin : Specify minimum depth(inclusive) of a sample to be included in analysis

When genotypes are filtered, they are marked as missing genotypes. 
Consequently, samples with missing genotype may or may not be included in the analysis.
That means samples with filtered genotypes may be dropped (`--impute drop`) 
or may still be included (`--impute mean` or `--impute hwe`). 
By default, genotypes are imputed to the mean value.
Please note that `--impute drop` is usually not recommended due to the incurred computational cost, 
as the null model may be estimated for each marker.

See next section about how you want to handle missing genotypes.

## Handle missing genotypes and phenotypes

When genotypes are missing (e.g. genotype = "./.") or genotypes are filtered out, 
there are three options to handle them: (1) impute to its mean(default option); (2) impute by HWE equilibrium; (3) remove from the model.
Use `--impute [mean|hwe|drop]` to specify which option to use.

When quantitative phenotypes are missing, for example, some samples have genotype files, but not phenotypes, 
rvtests can impute missing phenotype to its mean. 

*NOTE:* Do not use `--imputePheno` for binary trait.

In summary, the following two options can be used:

               --impute : Specify either of mean, hwe, and drop
          --imputePheno : Impute phenotype to mean by those have genotypes but no
                          phenotypes
                          
                          
## Specify groups (e.g burden unit)

Rare variants association tests are usually performed in groups of variants. 
The natural grouping unit is gene. Rvtests can read gene definition file in `refFlat` format,
and perform association for each gene. Use `--geneFile` option to specify the gene file name.
For example, `--geneFile refFlat_hg19.txt.gz` will use `refFlat_hg19.txt.gz` as gene definition file,
and then perform association tests for every gene. Use `--gene` to specify a subset of genes to test.
For example, `--gene CFH` will only test CFH gene.

Alternative grouping unit can be specified as *set*. 
These *sets* are treated similar to gene.
You can thus use `--setFile` to define sets (similar to `--geneFile` option), 
and use `--set` to define a specific set (similar to `--gene` option). 
Additionally, use `--setList` can specify a set to test from command line.

The format of a set file is: (1) set names; (2) ranges (e.g. chrom:begin-end);
For example, you have a set file, `example.set`, like this:

    set1 1:100-200,1:250-300
    set2 2:500-600
    
You can specify `--setFile example.set --set set2` to group variants 
within chromosome 2, position 500 to 600bp. 
If you want to test a particular region, for example, chromosome 2, position 500 to 550bp,
but do not want to make another file, you can use `--setList 2:500-600`.

In summary, options related to *Grouping Unit* are listed below: 

             --geneFile : specify a gene file (for burden tests)
                 --gene : specify which genes to test
              --setList : specify a list to test (for burden tests)
              --setFile : specify a list file (for burden tests, first two columns:
                          setName chr:beg-end)
                  --set : specify which set to test (1st column)


# Sex chromosome analysis

Rvtests support X chromosome analysis. In human X chromosome, there is PAR (pseudoautosomal region) and non-PAR region.
For males, there are two X allele in PAR region and one allele in non-PAR region.
While the PAR region is treated in the same way as autosomes, rvtests treat non-PAR region differently.
Below we will describe the details about how rvtests handles non-PAR region.

*Prepare data*. According to VCF standard, male genotype needs to coded as 0 or 1. For compatibility, rvtests also support 0/0 or 1/1 coding. 
In VCF files, male genotypes can be written as "0", "1", "0|0", "0/0", "1|1", "1/1". All other genotypes will be treated as missing.

*Genotype in regression models*. For consistence, male genotypes are converted to 0 or 2. When male dosages are provided, we expect the values in the VCF file are between 0.0 and 1.0, and then will model them as twice the value.

*MetaScore results*. If `--meta score` is specified, the output file `prefix.MetaScore.assoc.gz` includes both PAR-region and non-PAR region analysis. 
However, in the non-PAR region, the difference is that Hardy-Weinberg P-value and homozygous-reference/heterzygous/homozygous-alternative sample sizes are calculated **using female samples only**.

*Related individuals*. Just append `--xHemi` to the `vcf2kinship` (more details in [Kinship generation](#kinship-generation)) and `rvtest` command lines. Rvtests can recognize non-PAR region kinship file and use it in the analysis.

*PAR region*. PAR region is defined as two regions X:60001-2699520 and X:154931044-155260560. Use `--xLabel` can specify which chromosome has PAR region (default: 23 or X)
and use `--xParRegion` to specify PAR region (default: hg19, meaning '60001-2699520,154931044-155260560' in the UCSC build hg19, specify "hg18" will use PAR region definition in the UCSC build hg18, or specify "hg38" will use UCSC build 38).

# Kinship generation

Analysis of related individual usually requires estimation of kinship. You can a separate tool, `vcf2kinship`.
`vcf2kinship` is usually included in rvtests binary distribution or can be built from software source codes.

`vcf2kinship` can calculate pedigree kinship using a pedigree input file (PED format, see [Phenotype file](#phenotype-file), use option `--ped`).
The output file name is specified by `--prefix` option. If you use `--prefix output` then the output files will include `output.kinship`.

It can also calculate empirical kinship using genotype input file (VCF format, see [Genotype file (VCF)](#genotype-file-vcf), use option `--inVcf`).
For empirical kinship, you also need to specify the kinship model, either Balding-Nicols model (use option `--bn`) or Identity-by-state model (use option `--ibs`).

In sex chromosome analysis, it is often required to generate kinship on X chromosome regions, then you need to specify `--xHemi`. If your input VCF file has different X chromosome label (e.g. chromosome name is '23' instead of 'X'), you can use `--xLabel 23`.

If principal component decomposition (PCA) results are needed, you can use option `--pca`. Then output files with suffix '.pca' include PCA results.

When dealing with large input files, it is often preferred to use multiple CPU to speed up calculation using the option `--thread N` in which N is the number of CPU.

For example, to generate pedigree-based kinship (`--ped`) on both autosomal region and X chromosome (`--xHemi`) region, the command line is:

    vcf2kinship --ped input.ped --xHemi --out output

To generate empirical kinship (`--inVcf`) on both autosomal region and X chromosome (`--xHemi`) region using Balding-Nicols model, the command line is:

    vcf2kinship --inVcf input.vcf.gz --ped input.ped --bn --xHemi --out output

NOTE: you need to provide a pedigree file (PED) in the above case, as `vcf2kinship` needs the sex information of samples to construct kinship for sex chromosome.

For modern genetic datasets, genotype data is often stored by chromosomes, with each chromosome stored in a separate VCF file. In this case, kinship matrix can be first calcualted separately for each chromosome, and then combined using the python script [combineKinship.py](https://github.com/zhanxw/rvtests/tree/master/misc). Specifically, if you generated two kinship matrices for chr1 chr2 (chr1.kinship, and chr2.kinship), you can run the python script to combine them, i.e. 
```
python combineKinship.py -o prefix chr1.kinship chr2.kinship
```

The resulting kinship matrix is equivalent to the kinship matrix calculated using the merged vcf files.  

# Resources

## UCSC RefFlat Genes

refFlat_hg19.txt.gz

UCSC gene definition file in the refFlat format ([Details](http://genome.ucsc.edu/goldenpath/gbdDescriptionsOld.html#RefFlat)) in NCBI genome build 37 .

File link: <http://qbrc.swmed.edu/zhanxw/seqminer/data/refFlat_hg19.txt.gz>


## Gencode Genes

refFlat.gencode.v19.gz

Gencode gene definition version 19 in the refFlat format ([Details](http://www.gencodegenes.org/releases/19.html)). We have also previous versions of gene files and can provide upon request.

File link: <http://qbrc.swmed.edu/zhanxw/seqminer/data/refFlat.gencode.v19.gz>



# Frequently Asked Questions (FAQ)

* Does rvtests support binary traits of related-individuals?

Yes for meta-analysis models and no for most other models.
Proper analyses of related-individual are supported in meta-analysis model. You can refer to section [Meta-analysis tests](#meta-analysis-tests).
For many other association models, we notice that supporting binary traits for related individuals is complex and at current stage we have not found good solutions.

* Can you provide a list of command line options?

Rvtests have build-in help that can be found by executing `rvtest --help`.
We also put all available options in this [link](https://github.com/zhanxw/rvtests/wiki/Command-Line-Options).

* Can you provide standard error (SE) or confidence interval (CI) for the estimated Beta in the score model?

In the output of MetaScore model (--meta score), the standard error is the inverse of SQRT_V_STAT.
For example, if SQRT_V_STAT = 2, that means the standard error of estimated beta is 1/2 = 0.5.

* Why the INFORMATIVE_ALT_AC, N_REF and N_ALT columns have zero counts for certain chromosome X regions in meta-analysis models?

These counts are calculated from female individuals. If your study only has male samples, rvtests cannot report these counts. Because if a male carries a non-reference allele, we cannot conclude that this is heterozygous (0/1) site or homozygous alternatives (1/1) site.

* What happens when P-values equals to -1?

If rvtests fails to fit a certain model, it also fails to calculate P-value reliably. 
Rvtests will output P-value as -1 instead of any number between 0 and 1 to indicate that an error has occurred.
However, this should rarely happen. Please contact us if you have further questions.

* Why SKAT Q-values reported by rvtests are different from the SKAT R package?

We strictly follow the notations in the SKAT publication ([Wu et al. (2011) AJHG](http://www.hsph.harvard.edu/skat/)). However, in SKAT R package, its implementation is slightly different.
For example, in quantitative trait analysis, Q is divided by (2 * \hat{sigma2}) in the R package, but not in rvtests.
Although Q values can be different, the P-values from the two software packages should match (only in rare cases, numerical accuracy may cause minor differences).

* How rvtests handle multi-allelic variants?

In rvtests, we focus on bi-allelic variants, and thus by-default treat multi-allelic variants as bi-allelic variants. Any genotype that includes other than reference allele and the first alternative allele will be treated as missing.
For example, when the reference allele is 'A' in the VCF REF column and the alternative alleles are 'T,G' in the VCF ALT column, the genotype '0/2' will be treated as a missing genotype.

To properly analyze multi-allelic sites, we recommend coding multiple alternative alleles in separate lines. 
For example, a tri-allelic site A/T/G can be represented in two lines: 
the first line uses A as the reference allele, T as the alternative allele, and 0 or 1 to represent the frequency of alternative alleles; 
similary, the second line uses A as the reference allele, G as the alternative allele, and 0 or 1 to represent the frequency of alternative allele G.
The following table can be helpful to understand this coding scheme.

Genotype                     | A/A | A/T | A/G | T/T | T/G | G/G 
:----------------------------|:---:|:---:|:---:|:---:|:---:|:---
1st line (REF = A, ALT = T)  | 0/0 | 0/1 | 0/0 | 1/1 | 1/0 | 0/0 
2nd line (REF = A, ALT = G)  | 0/0 | 0/0 | 0/1 | 0/0 | 0/1 | 1/1 

The rationale of coding multi-allelic sites as multiple lines of bi-allelic sites is as follows:
(1) the majority of existing analysis software does not support multi-allelic sites coded in one line;
(2) imputation outputs (e.g. Michigan Imputation Server) multi-allelic sites in multiple lines.
If your genotype file cannot be encoded as multiple bi-allelic variants, rvtests has an experimental option `--multipleAllele`.
Specifying this option will enable the analysis of multi-allelic variants similar to the anlaysis of multiple bi-alleleic variants on-the-fly. 

* How to adjust for multiple comparisons?

Rvtests does not provide multiple comparison adjustment.
However, users can use p.adjust() function provided in R to calculate adjust p-values for multiple comparisons.
For example, users can choose "bonferroni" or "BH" (Benjamini and Hochberg).
For advanced users, qvalues can be calculated from the qvalue package provided by BioConductor.


# Feedback/Contact

Questions and requests can be sent to
Github issue page ([link](https://github.com/zhanxw/rvtests/issues))
or
Xiaowei Zhan ([zhanxw@umich.edu](mailto:zhanxw@umich.edu "mailto:zhanxw@umich.edu"))
or
Dajiang Liu ([dajiang.liu@outlook.com](mailto:dajiang.liu@outlook.com "mailto:dajiang.liu@outlook.com"))
or
Goncalo Abecasis ([goncalo@umich.edu](mailto:goncalo@umich.edu "mailto:goncalo@umich.edu")).


Rvtests is a collaborative effort from Xiaowei Zhan, Youna Hu, Bingshan Li, Dajiang Liu and Goncalo Abecasis.

We want to thank rvtests users and especially those who have provided valuable feedbacks. These users include: 
Xueling Sim, Scott Verize, Shuang Feng, Kevin Lu, Ruth Loos, Tessel Galesloot, Valerie Turcot, Stefan Gustafsson, Corbin Quick, Adam Locke, Michael Nalls, 
Jie Huang, Haitao Zhang, [the GIANT consortium](http://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium), [the GLGC consortium](http://lipidgenetics.org) and [the GSCAN consortium](https://ibg.colorado.edu/mediawiki/index.php/GSCAN).

