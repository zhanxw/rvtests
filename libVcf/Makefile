
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
   HOME = /home/zhanxw
endif
ifeq ($(UNAME), Darwin)
   HOME = /Users/zhanxw
endif

EXEC = lib-vcf.a
BASE = RangeList PeopleSet VCFUtil
OBJ = $(BASE:=.o)

GONCALO_INC = ../libsrc
GONCALO_LIB = ../libsrc/lib-goncalo.a

TABIX_INC = ../tabix
TABIX_LIB = ../tabix/libtabix.a

DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS -I. -I../base -I${GONCALO_INC} -I${TABIX_INC}

TARGET: CXXFLAGS = $(DEFAULT_CXXFLAGS)
TARGET: $(EXEC)

all: $(TARGET)

# 'make release' will: clean current build; build with -O4 flags;
# NOTE: don't use -j flag!
release: CXXFLAGS = -O4 $(DEFAULT_CXXFLAGS)
release: clean $(EXEC) 

debug: CXXFLAGS = -g $(DEFAULT_CXXFLAGS)
debug: $(EXEC)

$(EXEC): $(OBJ) 
	ar rcs $(EXEC) $(OBJ)

%: %.cpp # remove implicit Makefile rule of converting ABC.cpp to ABC
%.o: %.cpp %.h
	g++ -c $(CXXFLAGS) $<

clean:
	rm -rf $(EXEC) *.o
