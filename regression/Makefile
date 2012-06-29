LIB = lib-regression.a 
LIB-DBG = lib-dbg-regression.a 

SRC = LogisticRegression \
      LogisticRegressionScoreTest \
      LogisticRegressionPermutationTest \
      LinearRegression \
      LinearRegressionScoreTest \
      LinearRegressionPermutationTest \
      kbac_interface kbac \
      Skat qfc MixtureChiSquare \
      Table2by2

.PHONY: all release debug
all: release

GONCALO_INC = ../libsrc
# GONCALO_LIB = ../libsrc/lib-goncalo.a

EIGEN_INC = ../third/eigen
EIGEN_LIB =  # Eigen are header files only
$(EIGEN_INC) $(EIGEN_LIB):
	(cd ../third; make eigen)

GSL_INC = ../third/gsl/include

$(KBAC_LIB): 
	(cd kbac; make)

DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS #-Wall
release: CXXFLAGS = -O2 -DNDEBUG $(DEFAULT_CXXFLAGS)
release: $(LIB)
debug: CXXFLAGS = -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(LIB-DBG)


OBJ = $(SRC:%=%.o) $(KBAC_LIB)
$(LIB): $(OBJ)
	ar -cqs $@ $^
$(LIB-DBG): $(OBJ)
	ar -cqs $@ $^

define BUILD_each
  TAR := $(1).o
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): $$(SRC) 
	g++ -MMD $$(CXXFLAGS) -c $$< -I$(GONCALO_INC) -I$(EIGEN_INC) -I$(GSL_INC)
endef
$(foreach s, $(SRC), $(eval $(call BUILD_each, $(s))))


clean:
	rm -rf *.o *.a *.d

test: Test.cpp $(LIB)
	g++ -g $(CXXFLAGS) -c $< -I$(GONCALO_INC) 
	g++ -g $(CXXFLAGS) -o $@ Test.o $(LIB) $(GONCALO_LIB) -lz
	./test
