#!/usr/bin/env python
import sys, os, re

# convenient functions
def myopen(fn):
    import gzip
    f = gzip.open(fn)
    try:
        f.read(2)
        f.close()
        return gzip.open(fn)
    except:
        f.close()
        return open(fn)

def usage():
    print("%s -o prefix in1.kinship in2.kinship ..." % sys.argv[0] )

def getVariant(fn):
    ret = -1
    d = re.compile(r'Total \[ (\d+) \] variants are used to calculate autosomal kinship matrix.')
    for ln in myopen(fn):
        res = d.search(ln)
        if res:
            ret = res.groups()[0]
    return int(ret)

def loadKinship(fn):
    ids = []
    kin = []
    ncol = -1
    for i, ln in enumerate(myopen(fn)):
        fd = ln.strip().split()
        if i == 0:
            ncol = len(fd)
            continue
        ids.append(fd[:2])
        kin.append([float(i) for i in fd[2:]])
    nsample = ncol - 2
    if len(ids) != nsample or \
       any([nsample != len(i) for i in kin]):
        print >> sys.stderr, "Dimension not match in ", fn
    print >> sys.stderr, "Kinship file %s with %d samples loaded" % (fn, nsample)
    return ids, kin

def combineKinship(kin1, cnt1, kin2, cnt2):
    kin = kin1
    for i in xrange(len(kin1)):
        for j in xrange(len(kin1[0])):
            kin[i][j] = (kin1[i][j] * cnt1 + kin2[i][j] * cnt2) / (cnt1 + cnt2)
    cnt = cnt1 + cnt2
    print >> sys.stderr, "Combined %d variants from %d and %d variants" % (cnt, cnt1, cnt2)
    return kin, cnt

def outputKinship(ids, kin, prefix):
    outFn = open(prefix + '.kinship', 'w')
    # header
    header = ['FID', 'IID']
    header.extend([i[1] for i in ids])
    outFn.write('\t'.join(header))
    outFn.write('\n')
    # content
    nsample = len(ids)
    for i in xrange(nsample):
        content = [ids[i][0], ids[i][0]]
        content.extend(kin[i])
        outFn.write('\t'.join(map(str, content)))
        outFn.write('\n')
    outFn.close()
    
def isSameId(id1, id2):
    if len(id1) != len(id2):
        return False
    l = len(id1[0])
    for i in xrange(l):
        if id1[0][i] != id2[0][i]:
            return False
    l = len(id1[1])
    for i in xrange(l):
        if id1[1][i] != id2[1][i]:
            return False
    return True

if __name__ == '__main__':
    try:
        import getopt
        optlist, args = getopt.getopt(sys.argv[1:], 'o:')
        optlist = dict(optlist)
        fn = args
        prefix = optlist['-o']
    except:
        usage()
        raise
        sys.exit(1)

    # process log files
    logFiles = [i.replace('.kinship', '.vcf2kinship.log') for i in fn]
    numAuto = [getVariant(i) for i in logFiles]
    print >> sys.stderr, "Total %d variants from %d files are recognized." % (sum(numAuto), len(fn))
    
    # read first kinship
    ids, kin = loadKinship(fn[0])
    cnt = numAuto[0]
    
    # gradually combine rest kinships
    for i in xrange(1, len(fn)):
        ids2, kin2 = loadKinship(fn[i])
        cnt2 = numAuto[i]
        if not isSameId(ids, ids2):
            print >> sys.stderr, "Id not match, skippping"
            continue
        kin, cnt = combineKinship(kin, cnt, kin2, cnt2)

    # output results
    outputKinship(ids, kin, prefix = prefix)
    
