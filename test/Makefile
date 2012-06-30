check: check1 check2 check3 check4
check1:
	../executable/dbg/rvtest --inVcf test.vcf.gz --pheno try.test.ped --burden cmc --outputRaw --geneFile gene.txt --gene GENE1 --out check1
	diff -q check1.GENE1.6.3.data correct.check1.GENE1.6.3.data

check2:
	../executable/dbg/rvtest --inVcf test.vcf.gz --pheno try.test.ped --burden cmc --outputRaw --geneFile gene.txt --gene GENE1 --out check2 --impute drop
	diff -q check2.GENE1.5.3.data correct.check2.GENE1.5.3.data
check3: 
	../executable/dbg/rvtest --inVcf test.vcf.gz --pheno try.test.ped --burden cmc --outputRaw --geneFile gene.txt --gene GENE1 --out check3 --indvQualMin 3
	diff -q check3.GENE1.6.3.data correct.check3.GENE1.6.3.data
check4: 
	../executable/dbg/rvtest --inVcf test.vcf.gz --pheno try.test.ped --burden cmc --outputRaw --geneFile gene.txt --gene GENE1 --out check4 --indvDepthMin 3
	diff -q check4.GENE1.6.2.data correct.check4.GENE1.6.2.data

clean:
	rm -f check?.*
