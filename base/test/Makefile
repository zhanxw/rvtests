EXE = testLogger testIO testRangeList testUtils testRegex testSimpleMatrix testPedigree testKinship testTabixReader
all: $(EXE)
debug: all

CXX_FLAGS = -O0 -ggdb -I.. -I../../ -I../../base -I../../third/tabix/ \
-I../../third/pcre/include \
../../base/lib-dbg-base.a -lz -lbz2 \
../../third/pcre/lib/libpcreposix.a ../../third/pcre/lib/libpcre.a \
../../third/tabix/libtabix.a

define BUILD_each
  TAR := $(1)
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): $$(SRC) $(LIBS) 
	g++ -o $$@ $$< -MMD $$(CXX_FLAGS)
endef
$(foreach s, $(EXE), $(eval $(call BUILD_each, $(s))))

$(EXE): ../lib-dbg-base.a
check:
	./testIO && \
	./testLogger && \
	./testRangeList > testRangeList_output && \
	diff testRangeList_output testRangeList_output.correct && \
	./testUtils && \
	./testRegex && \
	./testSimpleMatrix  && \
	./testPedigree test.ped >& testPedigree.out && diff testPedigree.correct testPedigree.out && \
	./testKinship testKinship2.ped > testKinship2.out && diff testKinship2.out testKinship2.correct

#diff -q t1.output t1.correct

clean:
	-rm -f $(EXE) *.o *.d *.out

