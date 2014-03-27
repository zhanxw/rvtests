README
======

Here we present several exemplar command lines to use rvtests.
A complete manual can be found: 
[https://zhanxw.github.io/rvtests/](https://zhanxw.github.io/rvtests/)

= Run single variant analysis using Wald test

    ../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --out out1


= Specify the second column as the phenotype

    ../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --mpheno 2 --out out2


= Specify named phenotype y2

    ../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --pheno-name y2 --out out3


= Incorporate two named covariates: c1,c2

    ../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --covar covar --covar-name c1,c2 --out out4

= Handle covariate file with missing data (same syntax)

    ../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --covar covar.missing --covar-name c1,c2 --out out5

= Produce meta-analysis statistics from unrelated individuals (using covariates, regress covariates on pheontpye, then inverse normalize residuals, perform score test and generate covariance matrices)

    ../executable/rvtest --pheno pheno --inVcf example.vcf --meta score,cov --covar covar --covar-name c1,c2 --useResidualAsPhenotype --inverseNormal --out out6

= Generate kinship matrix from pedigree

    ../executable/vcf2kinship --ped pheno --bn --out output

= Generate empirical kinship matrix from VCF file using Balding-Nicols method

    ../executable/vcf2kinship --inVcf example.vcf --bn --out output

= Produce meta-analysis statistics from related individuals (using covariates, regress covariates on pheontpye, then inverse normalize residuals, perform score test and generate covariance matrices)

need to use a kinship matrix (see previous steps) output.kinship
then, use:

    ../executable/rvtest --pheno pheno --inVcf example.vcf --meta score,cov --covar covar --covar-name c1,c2 --useResidualAsPhenotype --inverseNormal --kinship output.kinship --out out7
n

= Produce meta-analysis statistics (like above), and also under dominant and recessive model

    ../executable/rvtest --pheno pheno --inVcf example.vcf --meta score,cov --covar covar --covar-name c1,c2 --useResidualAsPhenotype --inverseNormal --kinship output.kinship --out out8

= Burden tests (here shows CMC and Zeggini tests)

    ../executable/rvtest --pheno pheno --inVcf example.vcf.gz --setFile setFile --burden cmc,cmcWald,zeggini,zegginiWald --out out9

= Variable Threshold tests

    ../executable/rvtest --pheno pheno --inVcf example.vcf.gz --setFile setFile --vt price --out out10

= Kernel-based tests (SKAT and KBAC)

    ../executable/rvtest --pheno pheno --inVcf example.vcf.gz --setFile setFile --kernel kbac,skat --out out11
