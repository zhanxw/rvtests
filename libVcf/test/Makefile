SRC = t1 testVCFExtract testRangeList
all: $(SRC)
debug: all

CXX_FLAGS = -O0 -ggdb -I.. -I../../libsrc -I../../base -I../../third/tabix/ -I../../third/pcre/include ../../third/tabix/libtabix.a \
 ../../libVcf/lib-dbg-vcf.a \
 ../../base/lib-dbg-base.a ../../libsrc/lib-dbg-goncalo.a -lz -lbz2 \
../../third/pcre/lib/libpcreposix.a ../../third/pcre/lib/libpcre.a

t1: t1.cpp
	g++ -o $@ $< $(CXX_FLAGS)

testVCFExtract: testVCFExtract.cpp
	g++ -o $@ $< $(CXX_FLAGS)

testRangeList: testRangeList.cpp  ../../libVcf/lib-dbg-vcf.a 
	g++ -o $@ $< $(CXX_FLAGS)


check:
	./t1 > t1.output
	diff -q t1.output t1.correct

clean:
	-rm -f $(SRC)
