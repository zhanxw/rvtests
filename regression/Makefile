LIB = lib-regression.a 
LIB-DBG = lib-dbg-regression.a 

.PHONY: all release debug
all: release

GONCALO_INC = ../libsrc
# GONCALO_LIB = ../libsrc/lib-goncalo.a

EIGEN_INC = ../third/eigen
EIGEN_LIB =  # Eigen are header files only

# $(GONCALO_LIB):
# 	(cd ../libsrc; make)
# $(GONCALO_LIB_DBG):
# 	(cd ../libsrc; make debug)
$(EIGEN_INC) $(EIGEN_LIB):
	(cd ../third; make eigen)

DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS #-Wall
release: CXXFLAGS = -O2 $(DEFAULT_CXXFLAGS)
release: $(LIB)
debug: CXXFLAGS = -g -O0 $(DEFAULT_CXXFLAGS)
debug: $(LIB-DBG)

$(LIB): LogisticRegression.o LinearRegression.o Skat.o qfc.o Table2by2.o
	ar -cqs $@ $^
$(LIB-DBG): LogisticRegression.o LinearRegression.o Skat.o qfc.o Table2by2.o
	ar -cqs $@ $^

LogisticRegression.o: LogisticRegression.cpp LogisticRegression.h
	g++ $(CXXFLAGS) -c $< -I$(GONCALO_INC)
LinearRegression.o: LinearRegression.cpp LinearRegression.h
	g++ $(CXXFLAGS) -c $< -I$(GONCALO_INC)
Skat.o: Skat.cpp Skat.h
	g++ $(CXXFLAGS) -c $< -I$(GONCALO_INC) -I$(EIGEN_INC)
qfc.o: qfc.c
	g++ $(CXXFLAGS) -c $< -I$(GONCALO_INC)
Table2by2.o: Table2by2.cpp Table2by2.h
	g++ $(CXXFLAGS) -c $< -I$(GONCALO_INC)
clean:
	rm -rf *.o *.a *.d

test: Test.cpp $(LIB)
	g++ -g $(CXXFLAGS) -c $< -I$(GONCALO_INC) 
	g++ -g $(CXXFLAGS) -o $@ Test.o $(LIB) $(GONCALO_LIB) -lz
	./test
