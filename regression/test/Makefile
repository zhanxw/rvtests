SRC = testLinearRegression testLogisticRegression
all: $(SRC)

LIBS =  ../lib-dbg-regression.a ../../base/lib-dbg-base.a ../../libsrc/lib-dbg-goncalo.a
CXX_FLAGS = -O0 -ggdb -I.. -I../../libsrc -I../../base -I../../third/tabix/ -I../../third/pcre/include \
 $(LIBS) \
 ../../third/tabix/libtabix.a \
 ../../third/pcre/lib/libpcreposix.a ../../third/pcre/lib/libpcre.a \
 -lz -lbz2 

define BUILD_each
  TAR := $(1)
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): $$(SRC) $(LIBS) 
	g++ -o $$@ $$< -MMD $$(CXX_FLAGS)
endef
$(foreach s, $(SRC), $(eval $(call BUILD_each, $(s))))


#testLinearRegression: testLinearRegression.cpp
#	g++ -o $@ $< $(CXX_FLAGS)
#
#testLinearRegression: testLinearRegression.cpp
#	g++ -o $@ $< $(CXX_FLAGS)
#
check:
	./t1 > t1.output
	diff -q t1.output t1.correct

clean:
	-rm -f $(SRC)
