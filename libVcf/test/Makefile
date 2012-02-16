SRC = t1 testVCFExtract
all: $(SRC)


CXX_FLAGS = -O0 -ggdb -I.. -I../../libsrc -I../../base -I../../third/tabix/ -I../../third/pcre/include ../../third/tabix/libtabix.a \
 ../../base/lib-dbg-base.a ../../libsrc/lib-dbg-goncalo.a -lz -lbz2 \
../../third/pcre/lib/libpcreposix.a ../../third/pcre/lib/libpcre.a

t1: t1.cpp
	g++ -o $@ $< $(CXX_FLAGS)

testVCFExtract: testVCFExtract.cpp
	g++ -o $@ $< $(CXX_FLAGS)


check:
	./t1 > t1.output
	diff -q t1.output t1.correct

clean:
	-rm -f $(SRC)
