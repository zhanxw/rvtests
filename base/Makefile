UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
   HOME = /home/zhanxw
endif
ifeq ($(UNAME), Darwin)
   HOME = /Users/zhanxw
endif

LIB = lib-base.a
LIB_DBG = lib-dbg-base.a
BASE = Argument Exception IO OrderedMap Regex TypeConversion Utils Logger RangeList SimpleMatrix Pedigree Kinship
OBJ = $(BASE:=.o)

DEFAULT_CXXFLAGS = -I. -I../third/tabix -I../third/pcre/include

.PHONY: all release debug
all: release

# 'make release' will: clean current build; build with -O4 flags; backup executable files to release/ folder
# NOTE: don't use -j flag!
release: CXXFLAGS = -O4 $(DEFAULT_CXXFLAGS)
release: $(LIB)

debug: CXXFLAGS = -Wall -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(LIB_DBG)

$(LIB): $(OBJ)
	ar rcs $(LIB) $(OBJ)
$(LIB_DBG): $(OBJ)
	ar rcs $(LIB_DBG) $(OBJ)

%: %.cpp # remove implicit Makefile rule of converting ABC.cpp to ABC
%.o: %.cpp %.h
	g++ -c $(CXXFLAGS) $<

clean:
	rm -rf $(LIB) $(LIB_DBG) *.o *.d


# INC=-I. -I./cityhash/include -I/home/zhanxw/software/sparsehash/include
# LIB= -L./cityhash/lib -lcityhash

# map:
# 	g++ -std=c++0x -static -g -O0 -o test StringMap.cpp $(INC) $(LIB)
# mapTest:
# 	g++ -std=c++0x -static -O3 -o test StringMap.cpp $(INC) $(LIB) && xtime ./test
# mapTest2:
# 	g++ -std=c++0x -static -pg -g -o test StringMap.cpp $(INC) $(LIB) && xtime ./test
