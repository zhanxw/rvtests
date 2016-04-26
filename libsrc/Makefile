.PHONY: release debug
all: release

BASE = $(wildcard *.cpp)
OBJ = $(BASE:%.cpp=%.o)
OBJ_DBG = $(BASE:%.cpp=%_dbg.o)

release: lib-goncalo.a
lib-goncalo.a: $(OBJ)
	-rm $@
	ar -cqs lib-goncalo.a $^
%.o : %.cpp
	g++ -O4 -DNDEBUG -c $< -o $@

debug: lib-dbg-goncalo.a
lib-dbg-goncalo.a: $(OBJ_DBG)
	-rm $@
	ar -cqs lib-dbg-goncalo.a $^
%_dbg.o: %.cpp
	g++ -g -O0 -c $< -o $@

clean:
	-rm -rf *.o *.a
