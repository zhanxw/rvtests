for i in qt?.ped;
do
    python ped2vcf.py $i `basename $i .ped`
done;
