EXE = testLogger testIO testRangeList testUtils testRegex testSimpleMatrix \
      testPedigree testKinship testTabixReader testParRegion testKinshipToKinInbcoef \
      testCommonFunction testTypeConversion testSimpleTimer testProfiler testVersionChecker \
      testSocket testHttp
all: $(EXE)
debug: all

CXX_FLAGS = -O0 -ggdb -I.. -I../../ -I../../base -I../../third/gsl/include -I../../third/tabix/ \
-I../../third/pcre/include \
../../base/lib-dbg-base.a -lz -lbz2 \
../../third/pcre/lib/libpcreposix.a ../../third/pcre/lib/libpcre.a \
../../third/tabix/libtabix.a \
-DUSE_ACCURATE_TIMER \
-std=c++0x \
-Wall \
-lrt

define BUILD_each
  TAR := $(1)
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): $$(SRC) $(LIBS) 
	g++ -o $$@ $$< -MMD $$(CXX_FLAGS)
endef
$(foreach s, $(EXE), $(eval $(call BUILD_each, $(s))))

$(EXE): ../lib-dbg-base.a
check: kinship
	./testIO && \
	./testLogger && \
	./testRangeList > testRangeList_output && \
	diff testRangeList_output testRangeList_output.correct && \
	./testUtils && \
	./testRegex && \
	./testSimpleMatrix  && \
	./testPedigree test.ped >& testPedigree.out && diff testPedigree.correct testPedigree.out && \
        ./testParRegion && \
        ./testCommonFunction && \
        ./testTypeConversion > test.TypeConversion.output && diff test.TypeConversion.output test.TypeConversion.correct && \
        ./testSimpleTimer && \
        echo "All tests passed!"

kinship:
	./testKinship testKinship2.ped > testKinship2.out && diff testKinship2.out testKinship2.correct && \
        ./testKinshipToKinInbcoef a test.KinInbcoef.pedtest test.KinInbcoef.listtest > test.KinInbcoef.outtest && diff test.KinInbcoef.outtest test.KinInbcoef.outtest.correct && \
        ./testKinshipToKinInbcoef x test.KinInbcoefX.pedtestX2 test.KinInbcoefX.listtestX2 > test.KinInbcoefX.outtestX2 && diff test.KinInbcoefX.outtestX2 test.KinInbcoefX.outtestX2.correct

clean:
	-rm -f $(EXE) *.o *.d *.out

