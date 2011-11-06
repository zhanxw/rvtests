EXEC = rvtest

#STATGEN_LIB = ../statgen/lib/libStatGen.a ../statgen//lib/samtools-0.1.7a-hybrid/libbam.a
TABIX_LIB = ./tabix-0.2.2/libtabix.a

DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS

.PHONY: release debug
all: debug
release: CXXFLAGS = -O2 $(DEFAULT_CXXFLAGS)
release: $(EXEC)
debug: CXXFLAGS = -g $(DEFAULT_CXXFLAGS)
debug: $(EXEC)

rvtest: Main.cpp PeopleSet.h Utils.h RangeList.h OrderedMap.h IO.h Argument.h
	g++ -c $(CXXFLAGS) Main.cpp  -I. -D__ZLIB_AVAILABLE__ -lz -lbz2
	g++ -o $@ Main.o $(TABIX_LIB)  -lz -lbz2 -lm -lpcre -lpcreposix
clean: 
	rm -rf *.o $(EXEC)
doc: README
	pandoc README -o README.html

test: test1
test1: rvtest
	./rvtest --inVcf test.vcf --outVcf test1.out.vcf 
test2: rvtest
	./rvtest --inVcf test.vcf --outVcf test2.out.vcf --peopleIncludeID 1232,1455,1232 
test3: rvtest
	./rvtest --inVcf 100.vcf.gz --outVcf test3.vcf --peopleIncludeID 1160

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

# arg: Argument.h Argument.cpp
# 	g++ -g -o Argument Argument.cpp
# RangeList: RangeList_test.cpp RangeList.h RangeList_test_input
# 	g++ -c $(CXXFLAGS) RangeList_test.cpp -I../statgen/lib/include -I. -D__ZLIB_AVAILABLE__ -lz
# 	g++ -o $@ RangeList_test.o $(TABIX_LIB) $(STATGEN_LIB)  -lz -lm

# IO: IO_test.cpp IO.h 
# 	g++ -c $(CXXFLAGS) IO_test.cpp -I../statgen/lib/include -I. -D__ZLIB_AVAILABLE__ 
# 	g++ -o $@ IO_test.o $(TABIX_LIB) $(STATGEN_LIB)  -lz -lm -lbz2
