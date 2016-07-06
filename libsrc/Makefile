.PHONY: release debug
all: release

BASE = $(wildcard *.cpp)
OBJ = $(BASE:%.cpp=%.o)
OBJ_DBG = $(BASE:%.cpp=%_dbg.o)

release: CXX_FLAGS = -O4 -DNDEBUG
release: lib-goncalo.a
lib-goncalo.a: $(OBJ)
	-[ -f $@ ] && rm -f $@ || echo -n
	ar -cqs lib-goncalo.a $^
%.o : %.cpp
	g++ $(CXX_FLAGS) -c $< -o $@

debug: lib-dbg-goncalo.a
debug: CXX_FLAGS = -g -O0
lib-dbg-goncalo.a: $(OBJ_DBG)
	-[ -f $@ ] && rm -f $@ || echo -n
	ar -cqs lib-dbg-goncalo.a $^
%_dbg.o: %.cpp
	$(CXX) $(CXX_FLAGS) -c $< -o $@

profile: lib-dbg-goncalo.a
profile: CXX_FLAGS = -pg -O0

clean:
	-rm -rf *.o *.a
