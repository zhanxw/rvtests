all: release

ROOT=$(shell pwd)
export ROOT

include Makefile.common
include Makefile.lib

DIR_EXEC = ./executable
DIR_EXEC_DBG = ./executable/dbg

$(DIR_EXEC):
	mkdir -p $@
$(DIR_EXEC_DBG):
	mkdir -p $@

.PHONY: release debug lib lib-dbg clean tar doc

release: lib
	cd $(ROOT)/src && $(MAKE) release
	cd $(ROOT)/vcfUtils && $(MAKE) release

debug: lib-dbg
	cd $(ROOT)/src && $(MAKE) debug
	cd $(ROOT)/vcfUtils && $(MAKE) debug

##################################################
GitVersion.h: .git/HEAD .git/index
	-echo "const char *gitVersion = \"$(shell git rev-parse HEAD)\";" > $@
.git/HEAD .git/index:
	-echo "const char *gitVersion = \"not-a-git-repo\"" > GitVersion.h 

##################################################
## clean
##################################################
clean: 
	cd $(ROOT)/src && $(MAKE) clean
	cd $(ROOT)/vcfUtils && $(MAKE) clean

libclean:
	(cd $(ROOT)/base; make clean)
	(cd $(ROOT)/regression; make clean)
	(cd $(ROOT)/libVcf; make clean)

deepclean: clean libclean
	rm -rf *~
	(cd $(ROOT)/third; make deepclean)
	(cd $(ROOT)/libsrc; make clean)

# archive 
DATE=$(shell date '+%m%d')
tar:
	tar zvchf rvtest.$(DATE).tgz \
            Makefile .git/HEAD .git/index \
            *.h *.cpp \
            src third base libVcf regression libsrc

README.md: doc
doc: README.wiki
	java -jar third/wiki2html.jar README.wiki > README.html 
	pandoc -f html -t markdown README.html > README.md 
