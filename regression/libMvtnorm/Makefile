FC=gfortran

release: CXXFLAGS=-O2 -DNDEBUG
release: FCFLAGS=-O2 -DNDEBUG
debug: CXXFLAGS= -O0 -ggdb
debug: FCFLAGS=-O0 -ggdb

.SUFFIXES:

all: release
release: libMvtnorm.a
debug: libMvtnorm.a

libMvtnorm.a: mvtnorm.o mvt.o randomF77.o 
	ar rcs $@ $^ 

mvtnorm.o: mvtnorm.cpp mvtnorm.h
	$(CXX) $(CXXFLAGS) -c mvtnorm.cpp -o mvtnorm.o

randomF77.o: randomF77.c
	$(CC) $(CXXFLAGS) -c randomF77.c

mvt.o: mvt.f
	$(FC) $(FCFLAGS) -c mvt.f

clean:
	rm -f *.o *.a
