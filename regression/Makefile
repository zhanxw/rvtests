.PHONY: all release debug
all: release

ROOT?=..
include $(ROOT)/Makefile.common
include $(ROOT)/Makefile.lib

LIB = lib-regression.a 
LIB_DBG = lib-dbg-regression.a 

BASE = LogisticRegression \
       LogisticRegressionScoreTest \
       LogisticRegressionPermutationTest \
       LinearRegression \
       LinearRegressionScoreTest \
       LinearRegressionPermutationTest \
       LMM \
       kbac_interface kbac \
       Skat qfc MixtureChiSquare \
       Table2by2 \
       EigenMatrixInterface \
       FastLMM \
       GSLMinimizer \
       GSLIntegration \
       GrammarGamma \
       MetaCov \
       LogisticRegressionVT \
       LinearRegressionVT \
       FirthRegression \
       MultivariateNormalDistribution \
       MultivariateVT \
       FamSkat \
       cdflib \
       SkatO \
       MultipleTraitLinearRegressionScoreTest \
       FastMultipleTraitLinearRegressionScoreTest \
       Formula \
       ConjugateGradientSolver \
       BoltLMM

OBJ = $(BASE:%=%.o)
OBJ_DBG = $(BASE:%=%_dbg.o)

MVTNORM_OBJ = libMvtnorm/mvtnorm.o libMvtnorm/mvt.o libMvtnorm/randomF77.o
$(MVTNORM_OBJ):
	cd libMvtnorm; $(MAKE)

DEFAULT_CXXFLAGS ?= -D__STDC_LIMIT_MACROS
release: CXXFLAGS = -O2 -msse2 -DNDEBUG $(DEFAULT_CXXFLAGS)
release: $(LIB)
debug: CXXFLAGS = -Wall -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(LIB_DBG)
profile: CXXFLAGS = -Wall -pg -O0 $(DEFAULT_CXXFLAGS) -DEIGEN_NO_DEBUG  -DEIGEN_NO_STATIC_ASSERT 
profile: $(LIB_DBG)

$(LIB): $(OBJ) $(MVTNORM_OBJ)
	$(AR) -cqs $@ $^
$(LIB_DBG): $(OBJ_DBG) $(MVTNORM_OBJ)
	$(AR) -cqs $@ $^

define BUILD_each
  TAR := $(1).o
  TAR_DBG := $(1)_dbg.o
  SRC := $(1).cpp
  -include  $(1).d
  -include  $(1)_dbg.d
  $$(TAR) $$(TAR_DBG): $$(SRC) $$(EIGEN_INC) $$(GSL_INC) $$(NLOPT_INC)
	$$(CXX) -MMD $$(CXXFLAGS) -c $$< -o $$@ -I$(ROOT) -I$(GONCALO_INC) -I$(EIGEN_INC) -I$(GSL_INC) -I$(NLOPT_INC)
endef
$(foreach s, $(BASE), $(eval $(call BUILD_each, $(s))))

clean:
	rm -rf *.o *.a *.d

deepclean: clean
	(cd libMvtnorm; $(MAKE) clean)
