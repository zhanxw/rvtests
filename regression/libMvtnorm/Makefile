FC=gfortran

CXXFLAGS=-O2 -DNDEBUG
FCFLAGS=-O2 -DNDEBUG

.SUFFIXES:

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
