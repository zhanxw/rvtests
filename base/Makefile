.PHONY: all release debug
all: release

ROOT ?= ..
include $(ROOT)/Makefile.common
include $(ROOT)/Makefile.lib

LIB = lib-base.a
LIB_DBG = lib-dbg-base.a
BASE = Argument Exception IO OrderedMap Regex TypeConversion Utils Logger \
       RangeList SimpleMatrix Pedigree Kinship Profiler VersionChecker \
       Socket Http TextMatrix Indexer KinshipHolder RingMemoryPool
OBJ = $(BASE:%=%.o)
OBJ_DBG = $(BASE:%=%_dbg.o)

DEFAULT_CXXFLAGS ?= -std=c++0x -I. -I$(ROOT) $(addprefix -I, $(THIRD_INC))
# enable read over HTTP and FTP
DEFAULT_CXXFLAGS += -D_USE_KNETFILE

release: CXXFLAGS = -O2 $(DEFAULT_CXXFLAGS)
release: $(LIB)

debug: CXXFLAGS = -Wall -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(LIB_DBG)

profile: CXXFLAGS = -Wall -pg -O0 $(DEFAULT_CXXFLAGS)
profile: $(LIB_DBG)

$(LIB): $(OBJ)
	$(AR) rcs $@ $^
$(LIB_DBG): $(OBJ_DBG)
	$(AR) rcs $@ $^

define BUILD_each
  TAR := $(1).o
  TAR_DBG := $(1)_dbg.o
  SRC := $(1).cpp 
  -include  $(1).d
  -include  $(1)_dbg.d 
  $$(TAR) $$(TAR_DBG): $$(SRC) | $$(THIRD_INC)
	$$(CXX) -MMD $$(CXXFLAGS) $$(CXX_INCLUDE) -c $$< -o $$@
endef
$(foreach s, $(BASE), $(eval $(call BUILD_each, $(s))))

clean:
	rm -rf $(LIB) $(LIB_DBG) *.o *.d

