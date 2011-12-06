#!/usr/bin/python
import sys, os

def usage():
    print("%s in.ped outputPrefix" % sys.argv[0] )

def printVCFHeader(x):
    print "##VCF=4.0"
    print "#CHROM "+'\t'.join(x)
    
if __name__ == '__main__':
    inF, prefix = sys.argv[1:]
    data = [x.strip().split() for x in open(inF).xreadlines()]
    numPeople = len(data)
    numMarker = len(data[0]) - 6
    peopleName = [x[1] for x in data]
    printVCFHeader(peopleName)
    for i in xrange(6, numMarker):
        print "A"
