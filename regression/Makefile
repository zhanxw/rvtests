LIB = lib-regression.a 
GONCALO_INC = ../libsrc
GONCALO_LIB = ../libsrc/lib-goncalo.a

all: $(LIB)
$(LIB): LogisticRegression.o
	ar -cqs $(LIB) *.o
LogisticRegression.o: LogisticRegression.cpp LogisticRegression.h
	g++ -c $< -I$(GONCALO_INC)
clean:
	rm -rf *.o *.a *~
