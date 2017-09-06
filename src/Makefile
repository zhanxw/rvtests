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

BASE = Main \
      DataConsolidator \
      ModelFitter \
      TabixUtil \
      ModelManager \
      ModelParser \
      Model \
      GenotypeExtractor \
      VCFGenotypeExtractor \
      BGenGenotypeExtractor \
      KGGGenotypeExtractor \
      DataLoader \
      GenotypeCounter \

# Add the basename of your model files
# e.g. put SingleDummy on the right if your source codes are:
# SingleDummy.h and SingleDummy.cpp
BASE += SingleDummy

OBJ = $(BASE:%=%.o)
OBJ_DBG = $(BASE:%=%_dbg.o)

DEFAULT_CXXFLAGS ?= -std=c++0x -I. $(addprefix -I, $(THIRD_INC)) -DGIT_VERSION="\"$(GIT_VERSION)\"" -fopenmp
# enable read over HTTP and FTP
DEFAULT_CXXFLAGS += -D_USE_KNETFILE

release: CXX_FLAGS = -O2 -DNDEBUG $(DEFAULT_CXXFLAGS) $(STATIC_FLAG)
release: $(DIR_EXEC)/$(EXEC)
debug: CXX_FLAGS = -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(DIR_EXEC_DBG)/$(EXEC) 
profile: CXX_FLAGS = -pg -O0 $(DEFAULT_CXXFLAGS)
profile: $(DIR_EXEC_DBG)/$(EXEC) 

$(DIR_EXEC)/$(EXEC):  $(LIBS) \
                      $(OBJ) \
                      |$(DIR_EXEC)
	$(CXX) -o $@ $(OBJ) $(CXX_FLAGS) $(CXX_LIB)
$(DIR_EXEC_DBG)/$(EXEC): $(LIBS_DBG) \
                         $(OBJ_DBG) \
                         | $(DIR_EXEC_DBG)
	$(CXX) -o $@ $(OBJ_DBG) $(CXX_FLAGS) $(CXX_LIB_DBG)


##################################################
define BUILD_each
  TAR := $(1).o
  TAR_DBG := $(1)_dbg.o
  SRC := $(1).cpp 
  -include  $(1).d
  -include  $(1)_dbg.d 
  $$(TAR) $$(TAR_DBG): $$(SRC) | $$(THIRD_INC)
	$$(CXX) -MMD -c $$< -o $$@ $$(CXX_FLAGS) $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__
endef
$(foreach s, $(BASE), $(eval $(call BUILD_each, $(s))))

##################################################
## clean
##################################################
clean: 
	rm -rf *.o *.d \
        $(addprefix $(DIR_EXEC)/, $(EXEC)) \
        $(addprefix $(DIR_EXEC_DBG)/, $(EXEC))
libclean:
	(cd ../base; $(MAKE) clean)
	(cd ../regression; $(MAKE) clean)
	(cd ../libVcf; $(MAKE) clean)

deepclean: clean libclean
	rm -rf *~
	(cd ../third; $(MAKE) deepclean)
	(cd ../libsrc; $(MAKE) clean)
