all: release

ROOT?=..

include $(ROOT)/Makefile.common
include $(ROOT)/Makefile.lib

EXEC = rvtest

DIR_EXEC = $(ROOT)/executable
DIR_EXEC_DBG = $(ROOT)/executable/dbg

$(DIR_EXEC):
	mkdir -p $@
$(DIR_EXEC_DBG):
	mkdir -p $@

.PHONY: all release debug lib lib-dbg

release: CXX_FLAGS = -O2 -DNDEBUG $(DEFAULT_CXXFLAGS) $(STATIC_FLAG)
release: $(DIR_EXEC)/$(EXEC) 
$(DIR_EXEC)/$(EXEC):  $(LIBS) \
                      Main.o \
                      DataConsolidator.o \
                      ModelFitter.h \
                      TabixUtil.o \
                      |$(DIR_EXEC)
	$(CXX) -o $@ Main.o DataConsolidator.o TabixUtil.o $(CXX_FLAGS) $(CXX_LIB)

debug: CXX_FLAGS = -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(DIR_EXEC_DBG)/$(EXEC) 
$(DIR_EXEC_DBG)/$(EXEC): $(LIBS_DBG) \
                         Main.o \
                         DataConsolidator.o \
                         Collapsor.h \
                         ModelFitter.h \
                         TabixUtil.o \
                         | $(DIR_EXEC_DBG)
	$(CXX) -o $@ Main.o DataConsolidator.o TabixUtil.o $(CXX_FLAGS) $(CXX_LIB_DBG) 


##################################################
../GitVersion.h:
	(cd $(ROOT); make GitVersion.h)
-include Main.d
Main.o: Main.cpp ../GitVersion.h
	$(CXX) -MMD -c -o $@ $(CXX_FLAGS) $< $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__

-include DataConsolidator.d
DataConsolidator.o: DataConsolidator.cpp DataConsolidator.h
	$(CXX) -MMD -c -o $@ $(CXX_FLAGS) $< $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__

-include TabixUtil.d
TabixUtil.o: TabixUtil.cpp TabixUtil.h
	$(CXX) -MMD -c -o $@ $(CXX_FLAGS) $< $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__

##################################################
## clean
##################################################
clean: 
	rm -rf *.o *.d \
        $(addprefix $(DIR_EXEC)/, $(EXEC)) \
        $(addprefix $(DIR_EXEC_DBG)/, $(EXEC))
libclean:
	(cd base; make clean)
	(cd regression; make clean)
	(cd libVcf; make clean)
deepclean: clean libclean
	rm -rf *~
	(cd third; make deepclean)
	(cd libsrc; make clean)

##################################################
## tests
##################################################
test: test1
test1: $(DIR_EXEC_DBG)/$(EXEC)
	$(DIR_EXEC_DBG)/$(EXEC) --inVcf test.vcf.gz --outVcf test1.out.vcf 
test2: $(DIR_EXEC_DBG)/$(EXEC)
	$(DIR_EXEC_DBG)/$(EXEC) --inVcf test.vcf.gz --outVcf test2.out.vcf --peopleIncludeID P4,P2
test3: $(DIR_EXEC_DBG)/$(EXEC)
	$(DIR_EXEC_DBG)/$(EXEC) --inVcf test.vcf.gz --outVcf test3.vcf --peopleIncludeID P2,NotValid,P3 --peopleExcludeID P3
test4: $(DIR_EXEC_DBG)/$(EXEC)
	$(DIR_EXEC_DBG)/$(EXEC) --inVcf test.vcf.gz --make-bed test.plink

DajiangDataSet/qt1.vcf.gz: DajiangDataSet/qt1.ped
	(cd DajiangDataSet; bash cmd.sh);

DajiangDataSet := DajiangDataSet/qt1.vcf.gz DajiangDataSet/qt1.pheno

testSingle: rvtest $(DajiangDataSet)
	$(DIR_EXEC_DBG)/$(EXEC) --inVcf DajiangDataSet/qt1.vcf.gz --pheno DajiangDataSet/qt1.pheno --single score,wald
testBurden: rvtest $(DajiangDataSet)
	$(DIR_EXEC_DBG)/$(EXEC) --inVcf DajiangDataSet/qt1.vcf.gz --pheno DajiangDataSet/qt1.pheno --set DajiangDataSet/set.txt --burden cmc,zeggini,mb,exactCMC
testVt: rvtest $(DajiangDataSet)
	$(DIR_EXEC_DBG)/$(EXEC) --inVcf DajiangDataSet/qt1.vcf.gz --pheno DajiangDataSet/qt1.pheno --set DajiangDataSet/set.txt --vt cmc,zeggini,mb,skat
testKernel: rvtest $(DajiangDataSet)
	$(DIR_EXEC_DBG)/$(EXEC) --inVcf DajiangDataSet/qt1.vcf.gz --pheno DajiangDataSet/qt1.pheno --set DajiangDataSet/set.txt --kernel skat

# mem test:
testMemLeak: testMemLeak.cpp $(LIB)
	$(CXX) -g -O0 -o $@ $<  -I. $(INC)  $(LIB) -lz -lbz2 -lm -lpcre -lpcreposix

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

# # archive 
# DATE=$(shell date '+%m%d')
# tar:
# 	tar zvchf rvtest.$(DATE).tgz *.h *.cpp Makefile .git/HEAD .git/index third base libVcf regression libsrc

# arg: Argument.h Argument.cpp
# 	$(CXX) -g -o Argument Argument.cpp
# RangeList: RangeList_test.cpp RangeList.h RangeList_test_input
# 	$(CXX) -c $(CXXFLAGS) RangeList_test.cpp -I../statgen/lib/include -I. -D__ZLIB_AVAILABLE__ -lz
# 	$(CXX) -o $@ RangeList_test.o $(TABIX_LIB) $(STATGEN_LIB)  -lz -lm

# IO: IO_test.cpp IO.h 
# 	$(CXX) -c $(CXXFLAGS) IO_test.cpp -I../statgen/lib/include -I. -D__ZLIB_AVAILABLE__ 
# 	$(CXX) -o $@ IO_test.o $(TABIX_LIB) $(STATGEN_LIB)  -lz -lm -lbz2
# README.md:README.wiki
# doc: README.md
# 	java -jar third/wiki2html.jar README.wiki > README.html 
# 	pandoc -f html -t markdown README.html > README.md 
