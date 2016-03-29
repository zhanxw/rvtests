# UNAME := $(shell uname)
# ifeq ($(UNAME), Linux)
#    HOME = /home/zhanxw
# endif
# ifeq ($(UNAME), Darwin)
#    HOME = /Users/zhanxw
# endif

LIB = lib-base.a
LIB_DBG = lib-dbg-base.a
BASE = Argument Exception IO OrderedMap Regex TypeConversion Utils Logger \
       RangeList SimpleMatrix Pedigree Kinship Profiler VersionChecker \
       Socket Http
OBJ = $(BASE:%=%.o)

DEFAULT_CXXFLAGS = -std=c++0x -I. -I../third/tabix -I../third/pcre/include -I../third/bzip2
# enable read over HTTP and FTP
DEFAULT_CXXFLAGS += -D_USE_KNETFILE

.PHONY: all release debug
all: release

# 'make release' will: clean current build; build with -O4 flags; backup executable files to release/ folder
# NOTE: don't use -j flag!
release: CXXFLAGS = -O4 $(DEFAULT_CXXFLAGS)
release: $(LIB)

debug: CXXFLAGS = -Wall -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(LIB_DBG)

profile: CXXFLAGS = -Wall -pg -O0 $(DEFAULT_CXXFLAGS)
profile: $(LIB_DBG)

$(LIB): $(OBJ)
	$(AR) rcs $(LIB) $(OBJ)
$(LIB_DBG): $(OBJ)
	$(AR) rcs $(LIB_DBG) $(OBJ)

define BUILD_each
  TAR := $(1).o
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): $$(SRC) 
	$$(CXX) -MMD $$(CXXFLAGS) -c $$< $$(CXXFLAGS)
endef
$(foreach s, $(BASE), $(eval $(call BUILD_each, $(s))))

clean:
	rm -rf $(LIB) $(LIB_DBG) *.o *.d

