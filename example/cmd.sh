#!/bin/bash
# each command here is an example of how to use rvtests
# see README.md for more information

if [ ! -e ../executable/rvtest ]; then
    (cd ../; make)
fi

echo "Single variant analysis"
../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --out out1
../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --mpheno 2 --out out2
../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --pheno-name y2 --out out3
../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --covar covar --covar-name c1,c2 --out out4
../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --covar covar.missing --covar-name c1,c2 --out out5

echo "Meta-analysis (generate summary statistics)"
../executable/rvtest --pheno pheno --inVcf example.vcf --meta score,cov --covar covar --covar-name c1,c2 --useResidualAsPhenotype --inverseNormal --out out6
## for related individual, need to have their kinship estimated
../executable/vcf2kinship --ped pheno --bn --out output
../executable/vcf2kinship --inVcf example.vcf --bn --out output
## now use "--kinship" to specify kinship file
../executable/rvtest --pheno pheno --inVcf example.vcf --meta score,cov --covar covar --covar-name c1,c2 --useResidualAsPhenotype --inverseNormal --kinship output.kinship --out out7
../executable/rvtest --pheno pheno --inVcf example.vcf --meta score,cov,dominant,recessive --covar covar --covar-name c1,c2 --useResidualAsPhenotype --inverseNormal --kinship output.kinship --out out8
## binary related individuals are supported, binary phenotypes are automatically recognized (e.g. `y4` is a binary trait)
../executable/rvtest --pheno pheno --inVcf example.vcf --meta score,cov --kinship output.kinship --pheno-name y4 --out out7b
## NOT RECOMMENDED, this is just to demonstrate traits can be analyzed as quantitative traits using '--qtl'
../executable/rvtest --pheno pheno --inVcf example.vcf --meta score,cov --kinship output.kinship --pheno-name y4 --qtl --out out7q

echo "Rare-variant analysis"
../executable/rvtest --pheno pheno --inVcf example.vcf.gz --setFile setFile --burden cmc,cmcWald,zeggini,zegginiWald --out out9
../executable/rvtest --pheno pheno --inVcf example.vcf.gz --setFile setFile --vt price --out out10
../executable/rvtest --pheno pheno --inVcf example.vcf.gz --setFile setFile --kernel skat --out out11
../executable/rvtest --pheno pheno --inVcf example.vcf.gz --setFile setFile --kernel kbac,skat --pheno-name y4 --out out12

echo "Documentation of the examples can be found in README.md or https://github.com/zhanxw/rvtests/blob/master/README.md ."
