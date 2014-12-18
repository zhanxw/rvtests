LIB = lib-regression.a 
LIB-DBG = lib-dbg-regression.a 

BASE = LogisticRegression \
       LogisticRegressionScoreTest \
       LogisticRegressionPermutationTest \
       LinearRegression \
       LinearRegressionScoreTest \
       LinearRegressionPermutationTest \
       kbac_interface kbac \
       Skat qfc MixtureChiSquare \
       Table2by2 \
       EigenMatrixInterface \
       FastLMM \
       GSLMinimizer \
       GrammarGamma \
       MetaCov \
       LogisticRegressionVT \
       LinearRegressionVT \
       FirthRegression \
       MultivariateNormalDistribution \
       MultivariateVT \
       FamSkat \

OBJ = $(BASE:%=%.o)

.PHONY: all release debug
all: release

GONCALO_INC = ../libsrc
# GONCALO_LIB = ../libsrc/lib-goncalo.a

EIGEN_INC = ../third/eigen
EIGEN_LIB =  # Eigen are header files only
$(EIGEN_INC) $(EIGEN_LIB):
	(cd ../third; $(MAKE) eigen)

GSL_INC = ../third/gsl/include

MVTNORM_OBJ = libMvtnorm/mvtnorm.o libMvtnorm/mvt.o libMvtnorm/randomF77.o
$(MVTNORM_OBJ):
	cd libMvtnorm; $(MAKE)

DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS #-Wall
release: CXXFLAGS = -O2 -msse2 -DNDEBUG $(DEFAULT_CXXFLAGS)
release: $(LIB)
debug: CXXFLAGS = -Wall -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(LIB-DBG)
profile: CXXFLAGS = -Wall -pg -O0 $(DEFAULT_CXXFLAGS)
profile: $(LIB-DBG)

$(LIB): $(OBJ) $(MVTNORM_OBJ)
	$(AR) -cqs $@ $^
$(LIB-DBG): $(OBJ) $(MVTNORM_OBJ)
	$(AR) -cqs $@ $^

define BUILD_each
  TAR := $(1).o
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): $$(SRC) $$(EIGEN_INC)
	$$(CXX) -MMD $$(CXXFLAGS) -c $$< -I$(GONCALO_INC) -I$(EIGEN_INC) -I$(GSL_INC)
endef
$(foreach s, $(BASE), $(eval $(call BUILD_each, $(s))))

clean:
	rm -rf *.o *.a *.d

deepclean: clean
	(cd libMvtnorm; $(MAKE) clean)
