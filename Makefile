EXEC = rvtest

TABIX_INC = ./tabix-0.2.5
TABIX_LIB = ./tabix-0.2.5/libtabix.a

GONCALO_INC = ./libsrc
GONCALO_LIB = ./libsrc/lib-goncalo.a

VCF_INC = ./libVcf
VCF_LIB = ./libVcf/lib-vcf.a

REGRESSION_INC = ./regression
REGRESSION_LIB = ./regression/lib-regression.a

BASE_INC = ./base
BASE_LIB = ./base/lib-base.a

EIGEN_INC = ./eigen
EIGEN_LIB =  # Eigen are header files only

PCRE_INC = ./pcre/include
PCRE_LIB = ./pcre/lib/libpcre.a ./pcre/lib/libpcreposix.a

INC = -I$(TABIX_INC) -I$(REGRESSION_INC) -I$(VCF_INC) -I$(BASE_INC) -I$(EIGEN_INC) -I$(PCRE_INC) -I$(GONCALO_INC)
LIB = $(TABIX_LIB) $(REGRESSION_LIB) $(VCF_LIB) $(BASE_LIB) $(PCRE_LIB) $(GONCALO_LIB)

DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS #-Wall

.PHONY: release debug
all: debug
release: CXXFLAGS = -O2 $(DEFAULT_CXXFLAGS)
release: $(EXEC)
debug: CXXFLAGS = -g -O0 $(DEFAULT_CXXFLAGS)
debug: $(EXEC)

$(TABIX_LIB): tabix-0.2.5.tar.bz2 
	tar jvxf $<
	-(cd tabix-0.2.5; make;)

$(GONCALO_LIB): lib-goncalo.tgz 
	tar zvxf $<
	(cd libsrc; ./build.sh)

$(REGRESSION_LIB): $(shell ls -1 regression/*{cpp,h})
	(cd regression; make)

$(BASE_LIB): $(shell ls -1 base/*{cpp,h}) 
	(cd base; make)

$(VCF_LIB): $(shell ls -1 libVcf/*{cpp,h})
	(cd libVcf; make)

$(PCRE_LIB): pcre-8.21.tar.bz2 
	tar jvxf $<
	-(DIR=`pwd`; cd pcre-8.21; ./configure --prefix="$${DIR}"/pcre; make -j; make install)

rvtest: Main.o \
	VCFData.o \
	Collapsor.h ModelFitter.h \
	$(LIB)
	g++ -c $(CXXFLAGS) Main.cpp  -I. $(INC) -D__ZLIB_AVAILABLE__
	g++ -o $@ Main.o VCFData.o $(LIB) -lz -lbz2 -lm -lgsl -lblas

-include VCFData.d
VCFData.o: VCFData.cpp VCFData.h
	g++ -MMD -c $(CXXFLAGS) $<  -I. $(INC) -D__ZLIB_AVAILABLE__

-include Main.d
Main.o: Main.cpp
	g++ -MMD -c $(CXXFLAGS) $<  -I. $(INC) -D__ZLIB_AVAILABLE__

vcf2plink: vcf2plink.cpp $(LIB)
	#g++ -O4 -o $@ $<  -I. $(INC)  $(LIB) -lz -lbz2 -lm -lpcre -lpcreposix
	g++ -O0 -g -o $@ $<  -I. $(INC)  $(LIB) -lz -lbz2 -lm -lpcre -lpcreposix

plink2vcf: plink2vcf.cpp $(LIB)
	g++ -g -O0 -o $@ $<  -I. $(INC)  $(LIB) -lz -lbz2 -lm -lpcre -lpcreposix

vcf2merlin: vcf2merlin.cpp $(LIB)
	g++ -O4 -o $@ $<  -I. $(INC)  $(LIB) -lz -lbz2 -lm -lpcre -lpcreposix

vcf2ld: vcf2ld_window vcf2ld_neighbor
vcf2ld_window: vcf2ld_window.cpp $(LIB)
	g++ -O2 -g -o $@ $<  -I. $(INC)  $(LIB) -lz -lbz2 -lm -lpcre -lpcreposix

vcf2ld_neighbor: vcf2ld_neighbor.cpp $(LIB)
	g++ -O2 -g -o $@ $<  -I. $(INC)  $(LIB) -lz -lbz2 -lm -lpcre -lpcreposix

clean: 
	rm -rf *.o $(EXEC)

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
