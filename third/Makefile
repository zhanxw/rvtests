# download all third party library and compile them
.PHONY: tabix goncalo pcre eigen
all: tabix goncalo pcre eigen

tabix: tabix-0.2.5.tar.bz2
	-mv tabix-0.2.5.tar.bz2* tabix-0.2.5.tar.bz2
	tar jvxf tabix-0.2.5.tar.bz2
	ln -s -f tabix-0.2.5 tabix
	(cd tabix-0.2.5; make)

pcre: pcre-8.21.tar.gz
	-mv pcre-8.21.tar.gz* pcre-8.21.tar.gz
	tar zvxf pcre-8.21.tar.gz
	-(DIR=`pwd`; cd pcre-8.21; ./configure --prefix="$${DIR}"/pcre; make -j; make install)

#wget http://bitbucket.org/eigen/eigen/get/3.0.4.tar.bz2
eigen: eigen-3.0.4.tar.bz2
	tar jvxf eigen-3.0.4.tar.bz2
	-mv eigen-eigen-13a11181fc5a eigen-3.04
	ln -s -f eigen-3.04 eigen

clean:
	-cd tabix && make clean
#	-cd pcre && make clean
	-rm -f tabix eigen
