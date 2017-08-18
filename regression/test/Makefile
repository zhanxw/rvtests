EXE = testGSLIntegration \
      testLinearRegression \
      testLinearRegressionCov \
      testLogisticRegression \
      testLogisticRegressionCov \
      testLinearRegressionScoreTest \
      testLogisticRegressionScoreTest \
      testFirthRegression \
      testLogisticVT \
      testLinearVT \
      testMultivariateNormalDistribution \
      testLMM \
      testMixtureChiSquare \
      testMultipleTraitScoreTest \
      testFastMultipleTraitScoreTest \
      testFastMultipleTraitScoreTest.NA \
      testBoltLMM \
      testGSLMinimizer \
      testMatrixRef \

all: $(EXE)
debug: $(EXE)

../lib-dbg-regression.a: $(wildcard ../*.cpp) $(wildcard ../*.h)
	$(MAKE) -C .. debug

LIBS =  ../../libVcf/lib-dbg-vcf.a ../lib-dbg-regression.a ../../base/lib-dbg-base.a ../../libsrc/lib-dbg-goncalo.a
CXX_FLAGS = -fopenmp -O0 -ggdb -I../.. -I.. \
            -I../../libsrc -I../../base -I../../third/tabix/ \
            -I../../third/pcre/include \
            -I../../third/gsl/include \
            $(LIBS) \
            ../../third/tabix/libtabix.a \
            ../../third/pcre/lib/libpcreposix.a ../../third/pcre/lib/libpcre.a \
            ../../third/gsl/lib/libgsl.a \
            ../../third/nlopt/lib/libnlopt.a \
            -L../../third/zlib \
            -L../../third/bzip2 \
            -lz -lbz2 -lrt

define BUILD_each
  TAR := $(1)
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): $$(SRC) $(LIBS) 
	$$(CXX) -o $$@ $$< -MMD $$(CXX_FLAGS)
endef
$(foreach s, $(EXE), $(eval $(call BUILD_each, $(s))))

.PHONY: check
check: check1 check2 check3 check4 check5 check6 check8 check9
######################################################################
check1: output.R.lm output.cpp.lm
	python compare.py $^
output.R.lm: 
	Rscript testLinearRegression.R $@
output.cpp.lm: $(EXE) output.R.lm
	./testLinearRegression > $@
######################################################################
check2: output.R.lm.cov output.cpp.lm.cov
	python compare.py $^
output.R.lm.cov:
	Rscript testLinearRegressionCov.R $@
output.cpp.lm.cov: $(EXE) output.R.lm output.R.lm.cov
	./testLinearRegressionCov > $@
######################################################################
check3: output.R.glm output.cpp.glm
	python compare.py $^
output.R.glm:
	Rscript testLogisticRegression.R $@
output.cpp.glm: $(EXE)
	./testLogisticRegression > $@
######################################################################
check4: output.R.glm.cov output.cpp.glm.cov
	python compare.py $^
output.R.glm.cov:
	Rscript testLogisticRegressionCov.R $@
output.cpp.glm.cov: $(EXE)
	./testLogisticRegressionCov > $@
######################################################################
check5: output.score.R.lm output.score.cpp.lm
	python compare.py $^
output.score.R.lm:
	Rscript testLinearRegressionScoreTest.R $@
output.score.cpp.lm: $(EXE)
	./testLinearRegressionScoreTest > $@
######################################################################
check6: output.score.R.glm output.score.cpp.glm
	python compare.py $^
output.score.R.glm:
	Rscript testLogisticRegressionScoreTest.R $@
output.score.cpp.glm: $(EXE)
	./testLogisticRegressionScoreTest > $@

######################################################################
check7: correct.firth.out test.firth.out
	Rscript compare.R $^
correct.firth.out: 
	Rscript testFirthRegression.R $@
test.firth.out: $(EXE)
	./testFirthRegression > $@

######################################################################
check8: correct.logistic.mvt.out test.logistic.mvt.out
	Rscript compare.R $^
correct.logistic.mvt.out: 
	Rscript testLogisticVT.R $@
test.logistic.mvt.out: $(EXE)
	./testLogisticVT > $@

######################################################################
check9: correct.linear.mvt.out test.linear.mvt.out
	Rscript compare.R $^
correct.linear.mvt.out: 
	Rscript testLinearVT.R $@
test.linear.mvt.out: $(EXE)
	./testLinearVT > $@

######################################################################
check10: testMultipleTraitScoreTest.out \
         testFastMultipleTraitScoreTest.out \
         testFastMultipleTraitScoreTest.NA.out
	diff testMultipleTraitScoreTest.out correct.testMultipleTraitScoreTest.out
	diff testFastMultipleTraitScoreTest.out correct.testFastMultipleTraitScoreTest.out
	diff testFastMultipleTraitScoreTest.NA.out correct.testFastMultipleTraitScoreTest.NA.out
testMultipleTraitScoreTest.out : testMultipleTraitScoreTest
	./testMultipleTraitScoreTest > $@
testFastMultipleTraitScoreTest.out : testFastMultipleTraitScoreTest
	./testFastMultipleTraitScoreTest > $@
testFastMultipleTraitScoreTest.NA.out : testFastMultipleTraitScoreTest.NA
	./testFastMultipleTraitScoreTest.NA > $@

######################################################################
check11: testBoltLMM.out
	Rscript compare.bolt.corr.R $^ correct.testBoltLMM.out
testBoltLMM.out : $(EXE)
	./testBoltLMM $@

deepclean: clean
	-rm output.* input.*
clean:
	-rm -f $(EXE) *.o *.d
