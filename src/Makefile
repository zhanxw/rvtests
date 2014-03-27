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

release: CXX_FLAGS = -O2 -DNDEBUG $(DEFAULT_CXXFLAGS) $(STATIC_FLAG)
release: $(DIR_EXEC)/$(EXEC) 
$(DIR_EXEC)/$(EXEC):  $(LIBS) \
                      Main.o \
                      DataConsolidator.o \
                      ModelFitter.o \
                      TabixUtil.o \
                      |$(DIR_EXEC)
	$(CXX) -o $@ Main.o DataConsolidator.o ModelFitter.o TabixUtil.o $(CXX_FLAGS) $(CXX_LIB)

debug: CXX_FLAGS = -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(DIR_EXEC_DBG)/$(EXEC) 
profile: CXX_FLAGS = -pg -O0 $(DEFAULT_CXXFLAGS)
profile: $(DIR_EXEC_DBG)/$(EXEC) 
$(DIR_EXEC_DBG)/$(EXEC): $(LIBS_DBG) \
                         Main.o \
                         DataConsolidator.o \
                         Collapsor.h \
                         ModelFitter.o \
                         TabixUtil.o \
                         | $(DIR_EXEC_DBG)
	$(CXX) -o $@ Main.o DataConsolidator.o ModelFitter.o TabixUtil.o $(CXX_FLAGS) $(CXX_LIB_DBG) 


##################################################
../GitVersion.h:
	(cd $(ROOT); $(MAKE)  GitVersion.h)
-include Main.d
Main.o: Main.cpp ../GitVersion.h
	$(CXX) -MMD -c -o $@ $(CXX_FLAGS) $< $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__

-include DataConsolidator.d
DataConsolidator.o: DataConsolidator.cpp DataConsolidator.h
	$(CXX) -MMD -c -o $@ $(CXX_FLAGS) $< $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__

-include ModelFitter.d
ModelFitter.o: ModelFitter.cpp ModelFitter.h
	$(CXX) -MMD -c -o $@ $(CXX_FLAGS) $< $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__

-include TabixUtil.d
TabixUtil.o: TabixUtil.cpp TabixUtil.h
	$(CXX) -MMD -c -o $@ $(CXX_FLAGS) $< $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__

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

# # # archive 
# # DATE=$(shell date '+%m%d')
# # tar:
# # 	tar zvchf rvtest.$(DATE).tgz *.h *.cpp Makefile .git/HEAD .git/index third base libVcf regression libsrc

