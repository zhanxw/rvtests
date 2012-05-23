all: release
EXEC = rvtest
UTIL_EXEC = vcf2plink vcfSummary vcfConcordance # plink2vcf vcf2merlin vcf2ld_window vcf2ld_neighbor

DIR_EXEC = ./executable
DIR_EXEC_DBG = ./executable/dbg

$(DIR_EXEC):
	mkdir -p $@
$(DIR_EXEC_DBG):
	mkdir -p $@

##################################################
# Third-party libs.
TABIX_INC = third/tabix
TABIX_LIB = third/tabix/libtabix.a

EIGEN_INC = third/eigen
EIGEN_LIB =  # Eigen are header files only

PCRE_INC = third/pcre/include
PCRE_LIB = third/pcre/lib/libpcreposix.a third/pcre/lib/libpcre.a

$(TABIX_INC) $(TABIX_LIB):
	(cd third; make tabix)
$(EIGEN_INC) $(EIGEN_LIB):
	(cd third; make eigen)
$(PCRE_INC) $(PCRE_LIB):
	(cd third; make pcre)

THIRD_INC = $(TABIX_INC) $(EIGEN_INC) $(PCRE_INC)
THIRD_LIB = $(TABIX_LIB) $(PCRE_LIB)
##################################################
# Our libs.
BASE_INC = ./base
BASE_LIB = ./base/lib-base.a
BASE_LIB_DBG = ./base/lib-dbg-base.a

VCF_INC = ./libVcf
VCF_LIB = ./libVcf/lib-vcf.a
VCF_LIB_DBG = ./libVcf/lib-dbg-vcf.a

REGRESSION_INC = ./regression
REGRESSION_LIB = ./regression/lib-regression.a
REGRESSION_LIB_DBG = ./regression/lib-dbg-regression.a

GONCALO_INC = ./libsrc
GONCALO_LIB = ./libsrc/lib-goncalo.a
GONCALO_LIB_DBG = ./libsrc/lib-dbg-goncalo.a

$(BASE_LIB):
	(cd base; make)
$(BASE_LIB_DBG):
	(cd base; make debug)
$(VCF_LIB):
	(cd libVcf; make)
$(VCF_LIB_DBG):
	(cd libVcf; make debug)
$(REGRESSION_LIB): $(EIGEN_INC)
	(cd regression; make)
$(REGRESSION_LIB_DBG):
	(cd regression; make debug)
$(GONCALO_LIB):
	(cd libsrc; make)
$(GONCALO_LIB_DBG):
	(cd libsrc; make debug)

##################################################

INCLUDE = $(THIRD_INC) $(REGRESSION_INC) $(VCF_INC) $(BASE_INC) $(GONCALO_INC)
LIB = $(REGRESSION_LIB) $(VCF_LIB) $(BASE_LIB) $(GONCALO_LIB) $(THIRD_LIB) 
LIB_DBG = $(REGRESSION_LIB_DBG) $(VCF_LIB_DBG) $(BASE_LIB_DBG) $(GONCALO_LIB_DBG) $(THIRD_LIB)
CXX_INCLUDE = $(addprefix -I, $(INCLUDE))
CXX_LIB = $(LIB) -lz -lbz2 -lm -lgsl -lblas
CXX_LIB_DBG = $(LIB_DBG) -lz -lbz2 -lm -lgsl -lblas


DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS -std=c++0x #-Wall

.PHONY: release debug

lib: $(LIB)
lib-dbg: $(LIB_DBG)

release: CXX_FLAGS = -O2 $(DEFAULT_CXXFLAGS)
release: $(DIR_EXEC)/$(EXEC) util
$(DIR_EXEC)/$(EXEC): $(LIB) \
                     Main.o \
                     VCFData.o \
                     Collapsor.h ModelFitter.h \
                     |$(DIR_EXEC)
	g++ -o $@ Main.o VCFData.o $(CXX_FLAGS) $(CXX_LIB)

debug: CXX_FLAGS = -g -O0 $(DEFAULT_CXXFLAGS) 
debug: $(DIR_EXEC_DBG)/$(EXEC) util-dbg
$(DIR_EXEC_DBG)/$(EXEC): $(LIB_DBG) \
                         Main.o \
                         VCFData.o \
                         Collapsor.h ModelFitter.h \
                         | $(DIR_EXEC_DBG)
	g++ -o $@ Main.o VCFData.o $(CXX_FLAGS) $(CXX_LIB_DBG) 


##################################################
-include VCFData.d
VCFData.o: VCFData.cpp VCFData.h
	g++ -MMD -c $(CXXFLAGS) $< $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__

-include Main.d
Main.o: Main.cpp
	g++ -MMD -c $(CXXFLAGS) $< $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__


##################################################
# build utils
util: $(addprefix $(DIR_EXEC)/,$(UTIL_EXEC))
define BUILD_util
  TAR := $(DIR_EXEC)/$(notdir $(basename $(1)))
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): CXX_FLAGS = -O2 $(DEFAULT_CXXFLAGS)
  $$(TAR): $$(SRC) $(LIB) | $(DIR_EXEC)
	g++ -MMD -o $$@ $$< $$(CXX_FLAGS) $(CXX_INCLUDE) $(CXX_LIB)
endef
$(foreach s, $(UTIL_EXEC), $(eval $(call BUILD_util, $(s))))

util-dbg: $(addprefix $(DIR_EXEC_DBG)/,$(UTIL_EXEC))
define BUILD_util_dbg
  TAR := $(DIR_EXEC_DBG)/$(notdir $(basename $(1)))
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): CXX_FLAGS = -O0 -ggdb $(DEFAULT_CXXFLAGS)
  $$(TAR): $$(SRC) $(LIB_DBG) | $(DIR_EXEC_DBG)
	g++ -MMD -o $$@ $$< $$(CXX_FLAGS) $(CXX_INCLUDE) $(CXX_LIB_DBG)
endef
$(foreach s, $(UTIL_EXEC), $(eval $(call BUILD_util_dbg, $(s))))


clean: 
	rm -rf *.o *.d $(EXEC) \
        $(addprefix $(DIR_EXEC)/,$(UTIL_EXEC)) \
        $(addprefix $(DIR_EXEC_DBG)/,$(UTIL_EXEC))

deepclean: clean
	(cd base; make clean)
	(cd regression; make clean)
	(cd libVcf; make clean)

doc: README
	pandoc README -o README.html

test: test1
test1: rvtest
	./rvtest --inVcf test.vcf.gz --outVcf test1.out.vcf 
test2: rvtest
	./rvtest --inVcf test.vcf.gz --outVcf test2.out.vcf --peopleIncludeID P4,P2
test3: rvtest
	./rvtest --inVcf test.vcf.gz --outVcf test3.vcf --peopleIncludeID P2,NotValid,P3 --peopleExcludeID P3
test4: rvtest
	./rvtest --inVcf test.vcf.gz --make-bed test.plink

DajiangDataSet/qt1.vcf.gz: DajiangDataSet/qt1.ped
	(cd DajiangDataSet; bash cmd.sh);

DajiangDataSet := DajiangDataSet/qt1.vcf.gz DajiangDataSet/qt1.pheno

testSingle: rvtest $(DajiangDataSet)
	./rvtest --inVcf DajiangDataSet/qt1.vcf.gz --pheno DajiangDataSet/qt1.pheno --single score,wald
testBurden: rvtest $(DajiangDataSet)
	./rvtest --inVcf DajiangDataSet/qt1.vcf.gz --pheno DajiangDataSet/qt1.pheno --set DajiangDataSet/set.txt --burden cmc,zeggini,mb,exactCMC
testVt: rvtest $(DajiangDataSet)
	./rvtest --inVcf DajiangDataSet/qt1.vcf.gz --pheno DajiangDataSet/qt1.pheno --set DajiangDataSet/set.txt --vt cmc,zeggini,mb,skat
testKernel: rvtest $(DajiangDataSet)
	./rvtest --inVcf DajiangDataSet/qt1.vcf.gz --pheno DajiangDataSet/qt1.pheno --set DajiangDataSet/set.txt --kernel skat

# mem test:
testMemLeak: testMemLeak.cpp $(LIB)
	g++ -g -O0 -o $@ $<  -I. $(INC)  $(LIB) -lz -lbz2 -lm -lpcre -lpcreposix

# automated tests
autoTest: autoTest1 autoTest2
autoTest1: rvtest
	./rvtest --inVcf test.vcf.gz --outVcf test/try.test.vcf
	diff test/try.test.vcf test/correct.test.vcf

autoTest2: rvtest
	./rvtest --inVcf test.vcf.gz --make-bed test/try.test
	diff test/try.test.bim test/correct.test.bim
	diff test/try.test.fam test/correct.test.fam
	diff test/try.test.bed test/correct.test.bed

# archive 
DATE=$(shell date '+%m%d')
tar:
	tar zvchf rvtest.$(DATE).tgz *.h Main.cpp tabix*tar.bz2 

# arg: Argument.h Argument.cpp
# 	g++ -g -o Argument Argument.cpp
# RangeList: RangeList_test.cpp RangeList.h RangeList_test_input
# 	g++ -c $(CXXFLAGS) RangeList_test.cpp -I../statgen/lib/include -I. -D__ZLIB_AVAILABLE__ -lz
# 	g++ -o $@ RangeList_test.o $(TABIX_LIB) $(STATGEN_LIB)  -lz -lm

# IO: IO_test.cpp IO.h 
# 	g++ -c $(CXXFLAGS) IO_test.cpp -I../statgen/lib/include -I. -D__ZLIB_AVAILABLE__ 
# 	g++ -o $@ IO_test.o $(TABIX_LIB) $(STATGEN_LIB)  -lz -lm -lbz2
