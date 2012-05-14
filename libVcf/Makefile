UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
   HOME = /home/zhanxw
endif
ifeq ($(UNAME), Darwin)
   HOME = /Users/zhanxw
endif

.PHONY: all release debug
all: release 

LIB = lib-vcf.a
LIB_DBG = lib-dbg-vcf.a
BASE = RangeList PeopleSet VCFUtil VCFExtractor SiteSet
OBJ = $(BASE:=.o)

GONCALO_INC = ../libsrc
# GONCALO_LIB = ../libsrc/lib-goncalo.a

TABIX_INC = ../third/tabix
TABIX_LIB = ../third/tabix/libtabix.a

$(TABIX_INC) $(TABIX_LIB):
	(cd ../third; make tabix)

DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS -I. -I../base -I${GONCALO_INC} -I${TABIX_INC}

# 'make release' will: clean current build; build with -O4 flags;
# NOTE: don't use -j flag!
release: CXXFLAGS = -O4 $(DEFAULT_CXXFLAGS)
release: $(LIB) 

debug: CXXFLAGS = -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(LIB_DBG)

$(LIB): $(OBJ) 
	ar rcs $@ $^
$(LIB_DBG): $(OBJ) 
	ar rcs $@ $^

%: %.cpp # remove implicit Makefile rule of converting ABC.cpp to ABC
%.o: %.cpp %.h
	g++ -c $(CXXFLAGS) $<

clean:
	rm -rf $(LIB) $(LIB_DBG) *.o
