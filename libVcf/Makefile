.PHONY: all release debug
all: release 

LIB = lib-vcf.a
LIB_DBG = lib-dbg-vcf.a
BASE = PeopleSet VCFUtil PlinkInputFile PlinkOutputFile VCFInfo VCFInputFile VCFIndividual SiteSet VCFHeader # VCFExtractor
OBJ = $(BASE:=.o)

TABIX_INC = ../third/tabix
TABIX_LIB = ../third/tabix/libtabix.a

$(TABIX_INC) $(TABIX_LIB):
	(cd ../third; make tabix)

DEFAULT_CXXFLAGS = -std=c++0x -D__STDC_LIMIT_MACROS -I. -I../base -I${TABIX_INC}

# 'make release' will: clean current build; build with -O4 flags;
# NOTE: don't use -j flag!
release: CXXFLAGS = -O4 -DNDEBUG $(DEFAULT_CXXFLAGS)
release: $(LIB) 

debug: CXXFLAGS = -Wall -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(LIB_DBG)

$(LIB): $(OBJ) 
	ar rcs $@ $^
$(LIB_DBG): $(OBJ) 
	ar rcs $@ $^

%: %.cpp # remove implicit Makefile rule of converting ABC.cpp to ABC
%.o: %.cpp %.h
	g++ -c $(CXXFLAGS) $<

clean:
	rm -rf $(LIB) $(LIB_DBG) *.o *.d
