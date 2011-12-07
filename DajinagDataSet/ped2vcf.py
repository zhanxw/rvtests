#!/usr/bin/python
import sys, os

def usage():
    print("%s in.ped outputPrefix" % sys.argv[0] )

def printVCFHeader(x, f):
    f.write("##fileformat=VCFv4.0\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    f.write("\t")
    f.write('\t'.join(x))
    f.write("\n")

def geno2vcf(x):
    if x == '0':
	return '0/0'
    elif x == '1':
	return '0/1'
    elif x == '2':
	return '1/1'
    else:
	return './.'
    
if __name__ == '__main__':
    inF, prefix = sys.argv[1:]
    data = [x.strip().split() for x in open(inF).xreadlines()]
    numPeople = len(data)
    numMarker = len(data[0]) - 6 # PED file first 6 column: fid, pid, patid, matid, sex, pheno
    peopleName = [x[1] for x in data]

    # output VCF file
    fout = open(prefix + '.vcf', 'w')
    printVCFHeader(peopleName, fout)
    for idx, i in enumerate(xrange(numMarker)):
        fout.write('\t'.join(['1', str(idx+1), '.', 'A', 'T', '.', '.', '.', 'GT']))
	fout.write('\t')
	fout.write('\t'.join( (x[i+6] for x in data)))
	fout.write('\n')
    fout.close()

    # output phenotype file
    # col 5 (0-based) in .ped
    fout = open(prefix + '.pheno', 'w')
    for idx, i in enumerate(data):
	fout.write("%s\t%s\t%s\n" % (i[0], i[1], i[5]))
    fout.close()
    
		   
	
	
