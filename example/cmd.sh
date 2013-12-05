# each command here is an example of how to use rvtests
# see README for more information
../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --out out1
../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --mpheno 2 --out out2
../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --pheno-name y2 --out out3
../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --covar covar --covar-name c1,c2 --out out4
../executable/rvtest --pheno pheno --inVcf example.vcf --single wald --covar covar.missing --covar-name c1,c2 --out out5
../executable/rvtest --pheno pheno --inVcf example.vcf --meta score,cov --covar covar --covar-name c1,c2 --useResidualAsPhenotype --inverseNormal --out out6
../executable/vcf2kinship --ped pheno --bn --out output
../executable/vcf2kinship --inVcf example.vcf --bn --out output
../executable/rvtest --pheno pheno --inVcf example.vcf --meta score,cov --covar covar --covar-name c1,c2 --useResidualAsPhenotype --inverseNormal --kinship output.kinship --out out7
