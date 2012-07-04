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
            ../../third/gsl/lib/libgsl.a \
            -lz -lbz2 

define BUILD_each
  TAR := $(1)
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): $$(SRC) $(LIBS) 
	g++ -o $$@ $$< -MMD $$(CXX_FLAGS)
endef
$(foreach s, $(EXE), $(eval $(call BUILD_each, $(s))))

check: check1 check2 check3 check4
check1: $(EXE)
	Rscript testLinearRegression.R # output.R.lm
	./testLinearRegression > output.cpp.lm
	python compare.py output.R.lm output.cpp.lm
check2: $(EXE)
	Rscript testLinearRegressionCov.R # output.R.lm
	./testLinearRegressionCov > output.cpp.lm
	python compare.py output.R.lm output.cpp.lm
check3: $(EXE)
	Rscript testLogisticRegression.R # output.R.lm
	./testLogisticRegression > output.cpp.lm
	python compare.py output.R.lm output.cpp.lm
check4: $(EXE)
	Rscript testLogisticRegressionCov.R # output.R.lm
	./testLogisticRegressionCov > output.cpp.lm
	python compare.py output.R.lm output.cpp.lm

clean:
	-rm -f $(EXE) *.o *.d output.*
