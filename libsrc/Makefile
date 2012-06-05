.PHONY: release debug
all: release

SRC = *cpp *h

release: lib-goncalo.a
lib-goncalo.a: $(SRC)
	g++ -O4 -DNDEBUG -c *.cpp
	-rm $@
	ar -cqs lib-goncalo.a *.o
debug: lib-dbg-goncalo.a
lib-dbg-goncalo.a: $(SRC)
	g++ -g -O0 -c *.cpp
	-rm $@
	ar -cqs lib-dbg-goncalo.a *.o

clean:
	-rm -rf *.o *.a
