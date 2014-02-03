**Table of Contents**  

- [Introduction](#introduction)
- [Download](#download)
- [Quick Tutorial](#quick-tutorial)
	- [Single variant tests](#single-variant-tests)
	- [Groupwise tests](#groupwise-tests)
	- [Related individual tests](#related-individual-tests)
	- [Meta-analysis tests](#meta-analysis-tests)
- [Input files](#input-files)
	- [Genotype file (VCF)](#genotype-file-vcf)
	- [Phenotype file](#phenotype-file)
	- [Covariate file](#covariate-file)
	- [Trait transformation](#trait-transformation)
- [Models](#models)
	- [Single variant tests](#single-variant-tests-1)
	- [Burden tests](#burden-tests)
	- [Variable threshold models](#variable-threshold-models)
	- [Kernel models](#kernel-models)
	- [Utility models](#utility-models)
- [Association test options](#association-test-options)
	- [Sample inclusion/exclusion](#sample-inclusionexclusion)
	- [Variant site filters](#variant-site-filters)
	- [Genotyep filters](#genotyep-filters)
	- [Handle missing genotypes and phenotypes](#handle-missing-genotypes-and-phenotypes)
	- [Specify groups (e.g burden unit)](#specify-groups-eg-burden-unit)
- [Contact](#contact)

[![Build Status](https://travis-ci.org/zhanxw/rvtests.png?branch=master)](https://travis-ci.org/zhanxw/rvtests)

# Introduction

Rvtests, which stands for Rare Variant tests, is a flexible software package for genetic association studies. It is designed to support unrealted individual or related (family-based) individuals. Both quantitative trait and binary trait are supported. It includes a variety of association tests (e.g. single variant score test, burden test, variable threshold test, SKAT test, fast linear mixed model score test). It takes [VCF][vcf] format as genotype input file and takes PLINK format phenotype file and covariate file. From our practice, it is capable to analyze 8,000 related individuals using less than 400 Mb memory. 

[vcf]: http://www.1000genomes.com/

# Download

Source file can be downloaded from [github](https://github.com/zhanxw/rvtests/archive/master.zip) or [github page](https://github.com/zhanxw/rvtests).
Executable binary file can be downloaded from [here](http://www.sph.umich.edu/csg/zhanxw/software/rvtests/rvtests-latest.tar.gz).

# Quick Tutorial

Here is a quick example of how to use *rvtests* software in typical use cases.

## Single variant tests

    rvtests --inVcf input.vcf --pheno phenotype.ped --out output --single wald,score

This specifies single variant Wald and score test for association
tests for every variant in the `input.vcf` file. The 6th column of the phenotype file, `phenotype.ped`, which is in PLINK format, is used. Rvtests will automatically check whether the phenotype is binary trait or quantitative trait.
For binary trait, the recommended way of coding is to code controls as 1, cases as 2, missing phenotypes as -9 or 0.

For other types of association tests, you can refer to [Models](#Models)

## Groupwise tests
Groupwise tests includes three major kinds of tests.

* Burden tests: group variants, which are usually less than 1% or 5% rare variants, for association tests. The category includes: CMC test, Zeggini test, Madsen-Browning test, CMAT test, and rare-cover test.
* Variable threshold tests: group variants under different frequency thresholds.
* Kernel methods: suitable to tests rare variants having different directions of effects. These includes SKAT test and KBAC test. 

All above tests requires to group variants into a unit. The simplist case is to use gene as grouping unit. For different grouping method, see [Grouping](#Grouping). 

To perform rare variant tests by gene, you need to use `--geneFile` to specify the gene range in a refFlat format. We provided different gene definitions in the [Resources](#Resources) section. You can use `--gene` to specify which gene(s) to test. For example, specify `--gene CFH,ARMS2` will perform association tests on CFH and ARMS2 genes. If there is no providing `--gene` option, all genes will be tests.

The following command line demonstrate how to use CMC method, variable threshold method(proposed by Price) and kernel based method (SKAT by Shawn Lee and KBAC by
Dajiang Liu) to test every gene listed in *refFlat\_hg19\_uniq\_gene.txt.gz*.

    rvtests --inVcf input.vcf --pheno phenotype.ped --out output --geneFile refFlat_hg19_uniq_gene.txt.gz --burden cmc --vt price --kernel skat,kbac


## Related individual tests

To test related individuals, you will need to first create a kinship matrix:

    vcf2kinship --inVcf input.vcf --bn --out output

The option `--bn` means calculating empirical kinship using Balding-Nicols method. You can specifiy `--ibs` to obtain IBS kinship or use `--pedigree input.ped` to calculate kinship from known pedigree information.

Then you can use linear mixed model based association tests such as Fast-LMM score test, Fast-LMM LRT test and Grammar-gamma tests. An exemplar command is shown: 

    rvtests --inVcf input.vcf --pheno phenotype.ped --out output --kinship output.kinship --single famScore,famLRT,famGrammarGamma

## Meta-analysis tests

The meta-analysis models outputs association test results and genotype covariance matrix. These statistics can be used in rare variant association analysis.
We provide single variant score test and generate genotype covariance matrix. 
You can use command:
   
    rvtests --inVcf input.vcf --pheno phenotype.ped --meta score,cov --out output

In a more realistic scenario, you may want to adjust for covariates and want to inverse normalized residuals obtained in null model ([link](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.2852.html) to our methodology paper), then this command will work:
   
    rvtests --inVcf input.vcf --pheno phenotype.ped --covar example.covar --covar-name age,bmi --inverseNormal --useResidualAsPhenotype  --meta score,cov --out output


Here the `--covar` specify a covariate file, and `--covar-name` specify which covariates can used in the analysis. Covariate file format can be found [here](#Covariate file). `--inverseNormal --useResidualAsPhenotype` specifies trait transformation method. That means first fit a regression model of the phenotype on covariates (intercept automatically added), then the residuals are inverse normalized. Trait transformation details can be found [here](#Trait transformation).

We support both unrelated individuals and related indivudlas (e.g. family data). You need to append `--kinship input.kinship` to the command line:

    rvtests --inVcf input.vcf --pheno phenotype.ped --meta score,cov --out output --kinship input.kinship

The file `input.kinship` is calculated by `vcf2kinship` program, and usage to this program is described in [Related individual tests](#related-individual-tests).

Dominant and recessive disease models are supported by appending "dominant" and/or "recessive" after "--meta" option. For example, use "--meta dominant,recessive" will 
generate two files ".MetaDominant.assoc" and ".MetaRecessive.assoc". Details: In dominant models, genotypes 0/1/2 are coded as 0/1/1. In recessive models, genotypes 0/1/2 are 
coded as 0/0/1. Missing genotypes will be imputed to the mean.

# Input files

## Genotype file (VCF)

Rvtests supports VCF (Variant Call Format) files. Files in both plain txt format or gzipped format are supported. To use group-based rare variant tests, indexed the VCF files using [tabix](http://samtools.sourceforge.net/tabix.shtml) are required. 

Here are the commands to convert plain text format to bgzipped VCF format:

    (grep ^"#" $your_old_vcf; grep -v ^"#" $your_old_vcf | sed 's:^chr::ig' | sort -k1,1n -k2,2n) | bgzip -c > $your_vcf_file 
    tabix -f -p vcf $your_vcf_file

The above commands will (1) remove the `chr` prefix from chromosome names; (2) sort VCF files by chromosome first, then by chromosomal positions; (3) compress using bgzip; (4) create tabix index.

## Phenotype file

You can use `--mpheno $phenoypeColumnNumber` or `--pheno-name` to specify a given phenotype.

An example phenotype file, (`example.pheno`), has the following format: 

    fid iid fatid matid sex y1 y2 y3 y4
    P1 P1 0 0 0 1.7642934435605 -0.733862638327895 -0.980843608339726 2
    P2 P2 0 0 0 0.457111744989746 0.623297281416372 -2.24266162284447 1
    P3 P3 0 0 0 0.566689682543218 1.44136462889459 -1.6490100777089 1
    P4 P4 0 0 0 0.350528353203767 -1.79533911725537 -1.11916876241804 1
    P5 P5 0 0 1 2.72675074738545 -1.05487747371158 -0.33586430010589 2

Phenotype file is specified by the option `--pheno example.pheno` . The default phenotype column header is “`y1`”. If you want to use alternative columns as phenotype for association analysis (e.g the column with header y2), you may specific the header names using either

* --mpheno 2 
* --pheno-name y2

**NOTE:** to use “`--pheno-name`”, the  header line must starts with “`fid iid`” as PLINK requires.

In phenotype file, missing values can be denoted by NA or any non-numeric values. Individuals with missing phenotypes will be automatically dropped from subsequent association analysis. For each missing phenotype value, a warning will be generated and recorded in the log file.

## Covariate file

You can use `--covar` and `--covar-name` to specify covariates that will be used for single variant association analysis. This is an optional parameter. If you do not have covariate in the data, this option can be ignored. 

The covariate file, (e.g. `example.covar`) has a similar format as the phenotype file:

    fid iid fatid matid sex y1 y2 y3 y4
    P1 P1 0 0 0 1.911 -1.465 -0.817 1
    P2 P2 0 0 0 2.146 -2.451 -0.178 2
    P3 P3 0 0 0 1.086 -1.194 -0.899 1
    P4 P4 0 0 0 0.704 -1.052 -0.237 1
    P5 P5 0 0 1 2.512 -3.085 -2.579 1

The covariate file is specified by the `--covar` option (e.g. `--covar example.covar`). To specify covariates that will be used in the association analysis, the option `--covar-name` can be used. For example, when age, bmi and 3 PCs are used for association analysis, the following option can be specified for the rvtest program, i.e. 
`--covar example.covar --covar-name age,bmi,pc1,pc2,pc3`.

Note: Missing data in the covariate file can be labeled by any non-numeric value (e.g. NA). They will be automatically imputed to the mean value in the data file. 


## Trait transformation

In this meta-analysis, we use inversed normal transformed residuals in the association analysis, which is achieved by using a combination of `--inverseNormal`  and `--useResidualAsPhenotype`. Specifically, we first fit the null model by regressing phenotype on covariates. The residuals are then inverse normal transformed (see Appendix A more detailed formulae for transformation). Transformed residuals will be used to obtain score statistics. 

In meta analysis, an exemplar command for using rvtest looks like the following:

    ./rvtest --inVcf $vcf --pheno $example.pheno --covar example.covar --covar-name age,bmi --inverseNormal --useResidualAsPhenotype  --meta score,cov --out $output_prefix  

# Models
	
Rvtests support various association models.

## Single variant tests

Single variant | Model(*)    |Traits(#) | Covariates | Related / unrelated | Description
:--------------|:---------:|:------:|:----------:|:-------------------:|:-----------
Score test     |  score    |B, Q  |     Y      |         U           | Only null model is used to performed the test
Wald  test     |  wald     |B, Q  |     Y      |         U           | Only fit alternative model, and effect size will be estimated
Exact test     |  exact    |B     |     N      |         U           | Fisher's test
Fam LRT        |  famLRT   |Q     |     Y      |         R, U        | Fast-LMM model
Fam Score      |  famScore |Q     |     Y      |         R, U        | Fast-LMM model style likelihood ratio test
Grammar-gamma  |famGrammarGamma| Q     |     Y      |         R, U        | Grammar-gamma method


(*) Model columns list the regconized names in rvtests. For example, use `--single score` will apply score test.

(#) In trait column, B and Q stand for binary, quantitiave trait.


## Burden tests

Burden tests | Model(*)    |Traits(#) | Covariates | Related / unrelated | Description
:--------------|:---------:|:------:|:----------:|:-------------------:|:-----------
CMC             |  cmc       |B, Q  |     N      |         U           | Collapsing and combine rare variants by Bingshan Li.
Zeggini         |  zeggini   |B, Q  |     N      |         U           | Aggregate counts of rare variants by Morris Zeggini.
Madsen-Browning |  mb        |B     |     N      |         U           | Up-weight rare variant using inverse frequency from controls by Madsen.
Fp              |  fp        |B     |     N      |         U           | Up-weight rare variant using inverse frequency from controls by Danyu Lin.
Exact CMC       |  exactCMC  |B     |     N      |         U           | Collapsing and combine rare variants, then pefore Fisher's exact test.
RareCover       |  rarecover |B     |     N      |         U           | Find optimal grouping unit for rare variant tests by Thomas Hoffman.
CMAT            |  cmat      |B     |     N      |         U           | Test non-coding variants by Matt Z.
CMC Wald        |  cmcWald   |B, Q  |     N      |         U           | Collapsing and combine rare variants, then pefore Wald test.


(*) Model columns list the regconized names in rvtests. For example, use `--burden cmc` will apply CMC test.

(#) In trait column, B and Q stand for binary, quantitiave trait.


## Variable threshold models

Single variant | Model(*)    |Traits(#) | Covariates | Related / unrelated | Description
:--------------|:---------:|:------:|:----------:|:-------------------:|:-----------
Variable threshold model     |  vt    |B, Q  |     N      |         U           | Every rare-variant frequency cutoffs are tests by Alkes Price.  
Variable threshold CMC     |  cmc     |B, Q  |     N      |         U           | This models is natiive so that it output CMC test statistics under all possible frequency cutoffs.

(*) Model columns list the regconized names in rvtests. For example, use `--vt price` will apply score test.

(#) In trait column, B and Q stand for binary, quantitiave trait.



## Kernel models

Kernel | Model(*)    |Traits(#) | Covariates | Related / unrelated | Description
:--------------|:---------:|:------:|:----------:|:-------------------:|:-----------
SKAT     |  skat    |B, Q  |     Y      |         U           | Sequencing kernel association test by Shawn Lee.
KBAC     |  kbac     |B  |     N      |         U           | Kernel-based adaptive clustering model by Dajiang Liu.


(*) Model columns list the regconized names in rvtests. For example, use `--kernel skat` will apply SKAT test.
To further customize SKAT test, you can use *--kernel skat[nPerm=100:alpha=0.001:beta1=1:beta2=20]* to specify permutation counts, type-1 error, 
beta distribution parameters for upweighting rare variants. Rvtests will output a message showing: 

[INFO]  SKAT test significance will be evaluated using 10000 permutations at alpha = 0.001 (beta1 = 1.00, beta2 = 20.00)

(#) In trait column, B and Q stand for binary, quantitiave trait.


## Utility models


Rvtests has an usually option `--outputRaw`. When specify this, rvtests can output genotypes, phenotype, covariates(if any) and collapsed genotype to tabular files. These files can be imported into other software (e.g. R) for further analysis.


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
So to include sample A, B and C, you can provide a file, `people.txt`, looks like:

    A
    B
    C

Then use `--peopleIncludeFile people.txt` to include them in the analysis.


## Variant site filters

It is common that different frequency cutoffs are applied in rare-variant analysis.
Therefore, rvtests specify frequency cutoffs.

Frequency Cutoff

                 --freqUpper : Specify upper frequency bound to be included in analysis
                 --freqLower : Specify lower frequency bound to be included in analysis

Similar to sample inclusion/exclusion options, you can specify a range of variants to be included by 
specifying `--rangeList` option. For example `--rangeList 1:100-200` will include the chromosome 1 position 100bp to 200bp region.
Alternatively, use a separate file, `range.txt`, and `--rangeFile range.txt` to speicify association tests range.

                 --rangeList : Specify some ranges to use, please use chr:begin-end format.
                 --rangeFile : Specify the file containing ranges, please use chr:begin-end format.
                  --siteFile : Specify the file containing sites to include, please use "chr pos" format.

It is supported to filter variant site by site depth, minor allele count or annotation (annotated VCF file is needed).

              --siteDepthMin : Specify minimum depth(inclusive) to be incluced in analysis
              --siteDepthMax : Specify maximum depth(inclusive) to be incluced in analysis
                --siteMACMin : Specify minimum Minor Allele Count(inclusive) to be incluced in analysis
                  --annoType : Specify annotation type that is follwed by ANNO= in the VCF INFO field, regular expression is allowed

*NOTE*: `--annoType Nonsynonymous` will only analyze nonsynonymous variants where they have `ANNO=Nonsynonymous` in the INFO field. 
VCF with annotatino information are called annotated VCF here. And to annotate 
a VCF file, you can use [ANNO](https://github.com/zhanxw/anno), a fast and accurate annotation software.

## Genotyep filters

Genotype with low depth or low quality can be filtered out by these options:

              --indvDepthMin : Specify minimum depth(inclusive) of a sample to be incluced in analysis
              --indvDepthMax : Specify maximum depth(inclusive) of a sample to be incluced in analysis
               --indvQualMin : Specify minimum depth(inclusive) of a sample to be incluced in analysis

When genotypes are filtered, they are marked as missing genotypes. 
Consequently, samples with missing genotype may or may not be included in the analysis.
That means samples with genotypes may be dropped (`--impute drop`) 
or may still be included (`--impute mean` or `--impute hwe`). 
By default, genotypes are imputed to its means.
See next section about how you like to handle missing genotypes.


## Handle missing genotypes and phenotypes

When genotypes are missing (e.g. genotype = "./.") or gentoypes are filtered out, 
there are three options to handle them: (1) impute to its mean(default option); (2) impute by HWE equilibrium; (3) remove from the model.
Use `--impute [mean|hwe|drop]` to specify which option to use.

When quantitative phenotypes are missing, for example, some samples have gneotype files, but not phenotypes, 
rvtests can impute missing phenotype to its mean. 

*NOTE:* Do not use `--imputePheno` for binary trait.

In summary, the following two options can be used:

               --impute : Specify either of mean, hwe, and drop
          --imputePheno : Impute phenotype to mean by those have genotypes but no
                          phenotpyes
                          
                          
## Specify groups (e.g burden unit)

Rare variants association tests are usually performed in gruops of variants. 
The natural grouping unit is gene. Rvtests can read gene definition file in `refFlat` format,
and perform association for each gene. Use `--geneFile` option to specify the gene file name.
For example, `--geneFile refFlat_hg19.txt.gz` will use `refFlat_hg19.txt.gz` as gene definition file,
and then perform association tests for every gene. Use `--gene` to specify a subset of genes to test.
For example, `--gene CFH` will only test CFH gene.

Alternative grouping unit can be specified as *set*. 
These *sets* are treated similar to gene.
You can thus use `--setFile` to define sets (similar to `--geneFile` option), 
and use `--set` to define a specific set (similar to `--gene` option). 
Additionally, use `--setList` can speicify a set to test from command line.

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


# Contact

Questions and requests can be sent to Xiaowei Zhan
([zhanxw@umich.edu](mailto:zhanxw@umich.edu "mailto:zhanxw@umich.edu"))
or Goncalo Abecasis
([goncalo@umich.edu](mailto:goncalo@umich.edu "mailto:goncalo@umich.edu"))

Rvtests is a collaborative effort by Youna Hu, Bingshan Li, Dajiang Liu.





[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/zhanxw/rvtests/trend.png)](https://bitdeli.com/free "Bitdeli Badge")

