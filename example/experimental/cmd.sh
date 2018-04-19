#!/bin/bash
# This file includes example commands for experimental features
# Each command block shows one feature. Please note these feature may be changed in the future version.
# Please refer to README.md for all related information.

if [ ! -e ../../executable/rvtest ]; then
    (make -C ../..)
fi

echo "Support BGEN input genotype files"
echo "  replace --inVcf to --inBgen and --inBgenSample"
../../executable/rvtest --inBgen 1k.1m.bgen --inBgenSample 1k.1m.samples --pheno 1k.ped --out out1 --single score

echo "BGEN related utitlies"
echo "  Peek what is inside a BGEN file"
../../executable/bgenFileInfo --inBgen 1k.1m.bgen |head -n 42
echo "  Convert BGEN file to VCF file, with dosage, with hard coded genotypes and without variant id" # NOTE: BGEN has variant id and rs id per variant
../../executable/bgen2vcf --inBgen 1k.1m.bgen --inBgenSample 1k.1m.samples --showDS --hideVarId --out out.1k.1m

echo "Bolt-LMM model"
# NOTE: You need a binary PLINK file as inputs using --boltPlink
# Refer to ../README.md to understand how to prepare this binary PLINK file
../../executable/rvtest --pheno 1k.ped --inVcf 1k.1m.vcf.gz --meta bolt,boltCov --boltPlink 1k --out out2

echo "Documentation of the examples can be found in README.md or https://github.com/zhanxw/rvtests/blob/master/README.md ."
