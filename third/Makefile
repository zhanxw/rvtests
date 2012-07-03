# download all third party library and compile them
# .PHONY: tabix goncalo pcre eigen
all: tabix pcre eigen

tabix: tabix-0.2.5.tar.bz2
	-mv tabix-0.2.5.tar.bz2* tabix-0.2.5.tar.bz2
	tar jvxf tabix-0.2.5.tar.bz2
	ln -s -f tabix-0.2.5 tabix
	(cd tabix-0.2.5; make)

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


######################################################################
# download sources
eigen-3.1.0.tar.bz2:
	echo "obtain Eigen..."
	wget http://bitbucket.org/eigen/eigen/get/3.1.0.tar.bz2
	mv 3.1.0.tar.bz2 $@

gsl-1.15.tar.gz:
	echo "obtain GSL"
	wget http://mirrors.kernel.org/gnu/gsl/gsl-1.15.tar.gz

clean:
	-cd tabix && make clean
#	-cd pcre && make clean
	-rm -f tabix eigen
