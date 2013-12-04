all: release
release: util
debug: util-dbg

ROOT=..
include $(ROOT)/Makefile.common
include $(ROOT)/Makefile.lib

UTIL_EXEC = vcf2plink vcfSummary vcfConcordance \
            vcf2geno plink2vcf vcfExtractSite \
            vcf2ld_gene  \
            vcf2ld_window \
            vcfSummaryLite vcfAnnoSummaryLite vcfVariantSummaryLite vcfIndvSummary \
            vcfPair \
            vcf2kinship \
            vcfPeek \
            vcf2ld_neighbor 
#            gerp \
            # vcf2merlin 

DIR_EXEC = ../executable
DIR_EXEC_DBG = ../executable/dbg

$(DIR_EXEC):
	mkdir -p $@
$(DIR_EXEC_DBG):
	mkdir -p $@

.PHONY: release debug

# to build lib, we will use reverse order
# so third party lib is build first.
reverse = $(if $(1),$(call reverse,$(wordlist 2,$(words $(1)),$(1)))) $(firstword $(1))
lib: $(call reverse,$(LIB))
lib-dbg: $(call reverse,$(LIB_DBG))

##################################################
# build utils
##################################################
util: $(addprefix $(DIR_EXEC)/,$(UTIL_EXEC))
define BUILD_util
  TAR := $(DIR_EXEC)/$(notdir $(basename $(1)))
  SRC := $(addprefix ./, $(1).cpp)
  -include  $(1).d
  $$(TAR): CXX_FLAGS = -O2 -DNDEBUG $(DEFAULT_CXXFLAGS)  $(STATIC_FLAG)
  $$(TAR): $$(SRC) $(LIB) | $(DIR_EXEC)
	$(CXX) -MMD -o $$@ $$< $$(CXX_FLAGS) $(CXX_INCLUDE) $(CXX_LIB)
endef
$(foreach s, $(UTIL_EXEC), $(eval $(call BUILD_util, $(s))))

util-dbg: $(addprefix $(DIR_EXEC_DBG)/,$(UTIL_EXEC))
define BUILD_util_dbg
  TAR := $(DIR_EXEC_DBG)/$(notdir $(basename $(1)))
  SRC := $(addprefix ./, $(1).cpp)
  -include  $(1).d
  $$(TAR): CXX_FLAGS = -O0 -ggdb $(DEFAULT_CXXFLAGS)
  $$(TAR): $$(SRC) $(LIB_DBG) | $(DIR_EXEC_DBG)
	$(CXX) -MMD -o $$@ $$< $$(CXX_FLAGS) $(CXX_INCLUDE) $(CXX_LIB_DBG)
endef
$(foreach s, $(UTIL_EXEC), $(eval $(call BUILD_util_dbg, $(s))))


##################################################
## clean
##################################################
clean: 
	rm -rf *.o *.d \
        $(addprefix $(DIR_EXEC)/,$(UTIL_EXEC)) \
        $(addprefix $(DIR_EXEC_DBG)/,$(UTIL_EXEC))
