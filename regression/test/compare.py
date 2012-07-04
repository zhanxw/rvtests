#!/usr/bin/python
import sys, os

def usage():
    print("%s file1 file2: compare outputs" % sys.argv[0] )

if __name__ == '__main__':
    f1, f2 = sys.argv[1:3]
    a = [ln.strip().split() for ln in open(f1).xreadlines()]
    b = [ln.strip().split() for ln in open(f2).xreadlines()]
    if len(a) != len(b):
        print >> sys.stderr, "%s and %s differs in length" % (f1, f2)
        sys.exit(1)
    l = len(a)
    maxDiff = -999.9
    keyError = False
    for i in xrange(l):
        k1 = a[i][0]
        k2 = b[i][0]
        if k1 != k2:
            print >> sys.stderr, "%s and %s differs" % (k1, k2)
            keyError = True
            continue
        v1 = map(float, a[i][1:])
        v2 = map(float, b[i][1:])
        diff = max( [ abs(i-j) for i, j in zip(v1, v2) ])
        # print >> sys.stderr, "%s max diff is %f" % (k1, diff)
        if diff > maxDiff:
            maxDiff = diff
    if keyError:
        print "[WARNING] some results keys mismatch!"
    if maxDiff > 1e-2: # 10000
        print "[WARNING] maxDiff is HIGH: ", maxDiff
        sys.exit(1)
    else:
        print "[OK]"
