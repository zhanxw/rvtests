EXE = testReadVCF testReadVCFInfo testReadVCFByIndividual testReadVCFByRange \
      testReadMissing testPlinkInputFile testPlinkInputFile2 testVCFExtract \
      testTabixReader testBCFReader \
      testVCFMissingIndex testBCFMissingIndex \
      testReadBCFByRange \
      testVCFExtractChromX \
      testVCFSite \
      testVCFExtractChromXPar \
      testPlinkOutputFile \
      testPlinkOutputFile2 \
      testKGGInputFile

all: $(EXE)
debug: all

CXX_FLAGS = -O0 -ggdb -fopenmp \
  -I.. -I../.. -I../../libsrc -I../../base -I../../third/tabix/ -I../../third/pcre/include \
  -I../../third/samtools/bcftools \
  ../../libVcf/lib-dbg-vcf.a \
  ../../base/lib-dbg-base.a ../../libsrc/lib-dbg-goncalo.a \
  ../../third/tabix/libtabix.a \
  ../../third/pcre/lib/libpcreposix.a ../../third/pcre/lib/libpcre.a \
  -L ../../third/samtools/bcftools -lbcf\
  ../../third/samtools/libbam.a \
  -L ../../third/zlib \
  -L ../../third/bzip2 \
  -lz -lbz2 -lpthread

LIBS=../lib-dbg-vcf.a
lib:
	(cd ../../base; make debug)
	(cd ../; make debug)

define BUILD_each
  TAR := $(1)
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): $$(SRC) $(LIBS) 
	$$(CXX) -o $$@ $$< -MMD $$(CXX_FLAGS)
endef
$(foreach s, $(EXE), $(eval $(call BUILD_each, $(s))))

check: check1 check2 check3 check4 check5 check6 check7 check8 check9 check10 \
       check11 check12 check13 check14 check15 check16
check1:
	./testPlinkInputFile > testPlinkInputFile.output
	diff -q testPlinkInputFile.output testPlinkInputFile.output.correct
	./testPlinkInputFile2 > testPlinkInputFile2.output
	diff -q testPlinkInputFile2.output testPlinkInputFile2.output.correct
check2:
	./testVCFExtract > testVCFExtract.output
	diff -q testVCFExtract.output testVCFExtract.output.correct
check3:
	./testReadVCF > testReadVCF.output
	diff -q testReadVCF.output testReadVCF.output.correct
check4:
	./testReadVCFInfo > testReadVCFInfo.output
	diff -q testReadVCFInfo.output testReadVCFInfo.output.correct
check5:
	./testReadVCFByIndividual > testReadVCFByIndividual.output
	diff -q testReadVCFByIndividual.output testReadVCFByIndividual.output.correct
check6:
	./testReadVCFByRange > testReadVCFByRange.output
	diff -q testReadVCFByRange.output testReadVCFByRange.output.correct
check7:
	./testReadMissing > testReadMissing.output
	diff -q testReadMissing.output testReadMissing.output.correct
check8:
	./testTabixReader > testTabixReader.output
	diff -q testTabixReader.output testTabixReader.output.correct
check9:
	./testBCFReader > testBCFReader.output
	diff -q testBCFReader.output testBCFReader.output.correct
check10: 
	-./testVCFMissingIndex 2>&1 | tee >testVCFMissingIndex.output
	diff -q testVCFMissingIndex.output testVCFMissingIndex.output.correct
	-./testBCFMissingIndex 2>&1 | tee >testBCFMissingIndex.output
	diff -q testBCFMissingIndex.output testBCFMissingIndex.output.correct
check11:
	./testReadBCFByRange > testReadBCFByRange.output
	diff -q testReadBCFByRange.output testReadBCFByRange.output.correct
check12:
	./testVCFExtractChromX > testVCFExtractChromX.output
	diff -q testVCFExtractChromX.output testVCFExtractChromX.output.correct
check13:
	./testVCFSite > testVCFSite.output
	diff -q testVCFSite.output testVCFSite.correct
check14:
	./testVCFExtractChromXPar > testVCFExtractChromXPar.output
	diff -q testVCFExtractChromXPar.output testVCFExtractChromXPar.output.correct
check15:
	./testPlinkOutputFile > testPlinkOutputFile.output
	awk '{$$1 = NR; $$3 = 0; $$4 = 0; $$5 = 0; print ;}' testPlinkInputFile.fam | diff -b testPlinkOutputFile.output.fam -
	diff -q testPlinkOutputFile.output.bim testPlinkInputFile.bim
	diff -q testPlinkOutputFile.output.bed testPlinkInputFile.bed
	./testPlinkOutputFile2 > testPlinkOutputFile2.output
	-awk '{$$1 = NR; $$3 = 0; $$4 = 0; print ;}' testPlinkInputFile.fam | diff -b testPlinkOutputFile.output.fam -
	diff -q testPlinkOutputFile.output.bim testPlinkInputFile.bim
	diff -q testPlinkOutputFile.output.bed testPlinkInputFile.bed
check16:
	./testKGGInputFile > testKGGInputFile.vcf.output
	diff testKGGInputFile.vcf.output testKGGInputFile.vcf.output.correct
clean:
	-rm -f $(EXE) *.d *.output
