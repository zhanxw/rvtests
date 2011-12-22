#!/bin/bash
for i in qt?.ped;
do
    BASE=`basename $i .ped`
    python ped2vcf.py ${BASE}.ped ${BASE}
    ../tabix/bgzip -f ${BASE}.vcf
    ../tabix/tabix -f -p vcf ${BASE}.vcf.gz
done;
