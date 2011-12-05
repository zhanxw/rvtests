LIB = lib-regression.a 
GONCALO_INC = ../libsrc
GONCALO_LIB = ../libsrc/lib-goncalo.a

DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS -g #-Wall
release: CXXFLAGS = -O2 $(DEFAULT_CXXFLAGS)
release: $(LIB)
debug: CXXFLAGS = -g -O0 $(DEFAULT_CXXFLAGS)
debug: $(LIB)

all: $(LIB)
$(LIB): LogisticRegression.o LinearRegression.o
	ar -cqs $(LIB) *.o
LogisticRegression.o: LogisticRegression.cpp LogisticRegression.h
	g++ $(CXXFLAGS) -c $< -I$(GONCALO_INC)
LinearRegression.o: LinearRegression.cpp LinearRegression.h
	g++ $(CXXFLAGS) -c $< -I$(GONCALO_INC)
clean:
	rm -rf *.o *.a *~
test: Test.cpp $(LIB)
	g++ -g $(CXXFLAGS) -c $< -I$(GONCALO_INC) 
	g++ -g $(CXXFLAGS) -o $@ Test.o $(LIB) $(GONCALO_LIB) -lz
	./test
