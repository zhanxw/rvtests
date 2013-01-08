# download all third party library and compile them

LIB = tabix pcre eigen gsl bzip2
all: $(LIB)

tabix: tabix-0.2.6.tar.bz2
	-mv tabix-0.2.6.tar.bz2* tabix-0.2.6.tar.bz2
	tar jvxf tabix-0.2.6.tar.bz2
	ln -s -f tabix-0.2.6 tabix
	sed -i 's:-D_USE_KNETFILE::' tabix/Makefile
	sed -i 's:-D_USE_KNTFILE::' tabix/Makefile
	(cd tabix-0.2.6; make; make libtabix.so.1)

pcre: pcre-8.21.tar.gz
	-mv pcre-8.21.tar.gz* pcre-8.21.tar.gz
	tar zvxf pcre-8.21.tar.gz
	-(DIR=`pwd`; cd pcre-8.21; ./configure --prefix="$${DIR}"/pcre; make -j; make install)

gsl: gsl-1.15.tar.gz
	tar zvxf gsl-1.15.tar.gz
	-(DIR=`pwd`; cd gsl-1.15; ./configure --prefix="$${DIR}"/gsl; make -j; make install)

eigen: eigen-3.1.0.tar.bz2
	tar jvxf eigen-3.1.0.tar.bz2
	-mv eigen-eigen-ca142d0540d3 eigen-3.1.0
	ln -s -f eigen-3.1.0 eigen

bzip2: bzip2-1.0.6.tar.gz
	tar zvxf $<
	ln -s -f bzip2-1.0.6 bzip2
	(cd bzip2-1.0.6; make)

######################################################################
# download sources
eigen-3.1.0.tar.bz2:
	echo "obtain Eigen..."
	wget http://bitbucket.org/eigen/eigen/get/3.1.0.tar.bz2
	mv 3.1.0.tar.bz2 $@

gsl-1.15.tar.gz:
	echo "obtain GSL"
	wget http://mirrors.kernel.org/gnu/gsl/gsl-1.15.tar.gz

bzip2-1.0.6.tar.gz:
	echo "obtain Bzip2 library"
	wget http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
clean:
	-rm -rf $(LIB)
deepclean:
	-rm -rf $(LIB)
	-rm -rf tabix-0.2.6
	-rm -rf pcre-8.21
	-rm -rf gsl-1.15
	-rm -rf eigen-3.1.0
	-rm -rf bzip2-1.0.6

