EXE = t1 testVCFExtract testReadVCF
all: $(EXE)
debug: all

CXX_FLAGS = -O0 -ggdb -I.. -I../../libsrc -I../../base -I../../third/tabix/ -I../../third/pcre/include ../../third/tabix/libtabix.a \
 ../../base/lib-dbg-base.a ../../libsrc/lib-dbg-goncalo.a -lz -lbz2 \
../../third/pcre/lib/libpcreposix.a ../../third/pcre/lib/libpcre.a

define BUILD_each
  TAR := $(1)
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): $$(SRC) $(LIBS) 
	g++ -o $$@ $$< -MMD $$(CXX_FLAGS)
endef
$(foreach s, $(EXE), $(eval $(call BUILD_each, $(s))))

check:
	./t1 > t1.output
	diff -q t1.output t1.correct

clean:
	-rm -f $(EXE)
