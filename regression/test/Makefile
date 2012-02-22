EXE = testLinearRegression \
      testLinearRegressionCov \
      testLogisticRegression \
      testLogisticRegressionCov

all: $(EXE)
debug: $(EXE)

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
$(foreach s, $(EXE), $(eval $(call BUILD_each, $(s))))

check:
	./t1 > t1.output
	diff -q t1.output t1.correct

	Rscript testLinearRegression.R # output.R.lm
	./testLinearRegression > output.cpp.lm
	sleep 10

	Rscript testLogisticRegression.R # output.R.lm
	./testLogisticRegression > output.cpp.lm
	sleep 10

	Rscript testLinearRegressionCov.R # output.R.lm
	./testLinearRegressionCov > output.cpp.lm
	sleep 10

	Rscript testLogisticRegressionCov.R # output.R.lm
	./testLogisticRegressionCov > output.cpp.lm
	sleep 10

clean:
	-rm -f $(EXE) *.o *.d
