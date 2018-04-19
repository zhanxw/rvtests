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

.PHONY: release debug profile lib lib-dbg clean tar doc
release: lib
	$(MAKE) -C $(ROOT)/src release
	$(MAKE) -C $(ROOT)/vcfUtils release
	$(MAKE) -C $(ROOT)/bgenUtils release

debug: debug.rvt debug.vcfUtil debug.bgenUtil
debug.rvt: lib-dbg
	$(MAKE) -C $(ROOT)/src debug
debug.vcfUtil: lib-dbg
	$(MAKE) -C $(ROOT)/vcfUtils debug
debug.bgenUtil: lib-dbg
	$(MAKE) -C $(ROOT)/bgenUtils debug

profile: lib-dbg
	$(MAKE) -C $(ROOT)/src profile
	$(MAKE) -C $(ROOT)/vcfUtils profile
	$(MAKE) -C $(ROOT)/bgenUtils profile

##################################################
## clean
##################################################
clean: 
	$(MAKE) -C $(ROOT)/src clean
	$(MAKE) -C $(ROOT)/vcfUtils clean
	$(MAKE) -C $(ROOT)/bgenUtils clean

libclean:
	$(MAKE) -C $(ROOT)/base clean
	$(MAKE) -C $(ROOT)/regression clean
	$(MAKE) -C $(ROOT)/libVcf clean
	$(MAKE) -C $(ROOT)/libBgen clean

deepclean: clean libclean
	rm -rf *~
	$(MAKE) -C $(ROOT)/third deepclean
	$(MAKE) -C $(ROOT)/libsrc clean
	$(MAKE) -C $(ROOT)/base clean
	$(MAKE) -C $(ROOT)/regression deepclean
	$(MAKE) -C $(ROOT)/libVcf clean
	$(MAKE) -C $(ROOT)/libBgen clean

# archive 
DATE=$(shell date '+%m%d')
tar:
	tar zvchf rvtest.$(DATE).tgz \
            Makefile .git/HEAD .git/index \
            *.h *.cpp \
            src third base libVcf libBgen regression libsrc

wiki2md: README.wiki
	java -jar third/wiki2html.jar README.wiki > README.html 
	pandoc -f html -t markdown README.html > README.md 

md2wiki: README.md
	pandoc -f markdown_github+pipe_tables -t mediawiki README.md > README.wiki
