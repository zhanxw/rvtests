EXE = testBGenFile \
      testBGenIndex \
      testBGenVariant \
      testBGenFileByRange

all: $(EXE)
debug: all

CXX_FLAGS = -O0 -ggdb -fopenmp \
  -I.. -I../.. -I../../libsrc -I../../base -I../../third/tabix/ -I../../third/pcre/include \
  ../../libBgen/lib-dbg-bgen.a \
  ../../base/lib-dbg-base.a ../../libsrc/lib-dbg-goncalo.a \
  ../../third/pcre/lib/libpcreposix.a ../../third/pcre/lib/libpcre.a \
  ../../third/zstd/lib/libzstd.a ../../third/sqlite/lib/libsqlite3.a \
  ../../third/samtools/libbam.a \
  ../../third/zlib/libz.a \
  ../../third/bzip2/libbz2.a \
  -lpthread -ldl

LIBS=../lib-dbg-bgen.a
lib: $(LIBS)
$(LIBS):
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

check: check1 check2 check3 
check1: testBGenFile
	# BGEN v1.1
	./testBGenFile example.v11.bgen > example.v11.bgen.vcf && diff -q example.v11.bgen.vcf example.v11.bgen.vcf.correct
	# BGEN complex alleles
	./testBGenFile complex.bgen > complex.bgen.vcf && diff -q complex.bgen.vcf complex.bgen.vcf.correct
	# BGEN variable bits representation
	./testBGenFile complex.1bits.bgen >complex.1bits.bgen.vcf && diff -q complex.1bits.bgen.vcf complex.1bits.bgen.vcf.correct
	./testBGenFile complex.31bits.bgen >complex.31bits.bgen.vcf && diff -q complex.31bits.bgen.vcf complex.1bits.bgen.vcf.correct
	# BGEN ZSTD compress
	./testBGenFile example.16bits.zstd.bgen > example.16bits.zstd.bgen.vcf && diff -q example.16bits.zstd.bgen.vcf example.16bits.zstd.bgen.vcf.correct
check2: testBGenIndex
	./testBGenIndex example.16bits.bgen.bgi > testBgenIndex.output
	diff -q testBgenIndex.output testBgenIndex.output.correct
check3: testBGenVariant
	./testBGenVariant > testBGenVariant.output
	diff -q testBGenVariant.output testBGenVariant.output.correct
check4: testBGenFileByRange
	#./testBGenFile example.16bits.bgen |awk '$0 ~ /^#/ || ($2  >=2000 && $2 < 5000) || ($2 >=6000 && $2 < 8000)' |sort -k1,1r -k2,2g > testBGenFileByRange.output.correct
	./testBGenFileByRange example.16bits.bgen > testBGenFileByRange.output
	diff -q testBGenFileByRange.output testBGenFileByRange.output.correct

clean:
	-rm -f $(EXE) *.d *.output *.vcf
