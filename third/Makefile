# download all third party library and compile them

LIB = pcre eigen gsl bzip2 zlib rapidjson sqlite zstd nlopt libdeflate htslib
all: $(LIB)

pcre: pcre-8.36.tar.gz
	-mv pcre-8.36.tar.gz* pcre-8.36.tar.gz
	tar zxf pcre-8.36.tar.gz
	-(DIR=`pwd`; cd pcre-8.36; ./configure --prefix="$${DIR}"/pcre; make -j; make install)

gsl: gsl-1.16.tar.gz
	tar zxf gsl-1.16.tar.gz
	-(DIR=`pwd`; cd gsl-1.16; ./configure --prefix="$${DIR}"/gsl; make -j; make install)

eigen: eigen-3.3.7.tar.bz2
	-rm -rf eigen-eigen*
	tar jxf eigen-3.3.7.tar.bz2
	-mv eigen-eigen-* eigen-3.3.7
	ln -s -f eigen-3.3.7 eigen

zlib: zlib-1.2.8.tar.gz
	tar zxf $<
	ln -s -f zlib-1.2.8 zlib
	(cd zlib; ./configure; make)

bzip2: bzip2-1.0.6.tar.gz
	tar zxf $<
	ln -s -f bzip2-1.0.6 bzip2
	(cd bzip2-1.0.6; make)

nlopt: nlopt-2.4.2.tar.gz
	tar zxf $<
	-(DIR=`pwd`; cd nlopt-2.4.2; ./configure --prefix="$${DIR}"/nlopt; make -j; make install)

rapidjson: rapidjson-1.0.2.tar.gz
	tar zxf $<
	ln -s -f rapidjson-1.0.2 rapidjson
	sh -c '# no need to compile, just need its header'

sqlite: sqlite-autoconf-3200100.tar.gz
	tar zxf $<
	-(DIR=`pwd`; cd sqlite-autoconf-3200100; ./configure --prefix="$${DIR}"/sqlite --enable-dynamic-extensions=no; make -j; make install)

zstd: zstd-1.3.1.tar.gz
	tar zxf $<
	ln -s -f zstd-1.3.1 zstd
	(cd zstd; make)

libdeflate: libdeflate-1.3.tar.gz
	tar zxf $<
	(DIR=`pwd`; cd libdeflate-1.3; make install PREFIX="$${DIR}"/libdeflate)

htslib: htslib-1.9.tar.gz libdeflate
	tar zxf $<
	(DIR=`pwd`; cd htslib-1.9; autoconf ; ./configure --prefix="$${DIR}"/htslib --with-libdeflate CPPFLAGS=-I"$${DIR}"/libdeflate/include LDFLAGS=-L"$${DIR}"/libdeflate/lib --disable-libcurl; make install)

######################################################################
# download sources
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
  DOWNLOAD = wget --no-check-certificate $(1) 
  DOWNLOAD_RENAME = wget --no-check-certificate $(1) -O $(2)
else
  DOWNLOAD = curl -L $(1) -O
  DOWNLOAD_RENAME = curl -L $(1) -o $(2)
endif

pcre-8.36.tar.gz:
	echo "obtain PCRE..."
	$(call DOWNLOAD,ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.36.tar.gz)

eigen-3.3.7.tar.bz2:
	echo "obtain Eigen..."
	$(call DOWNLOAD_RENAME,http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2,$@)

gsl-1.16.tar.gz:
	echo "obtain GSL"
	$(call DOWNLOAD,http://mirrors.kernel.org/gnu/gsl/gsl-1.16.tar.gz)

zlib-1.2.8.tar.gz:
	echo "obtain zlib library"
	$(call DOWNLOAD,http://zlib.net/zlib-1.2.8.tar.gz)

bzip2-1.0.6.tar.gz:
	echo "obtain Bzip2 library"
	$(call DOWNLOAD,http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz)

nlopt-2.4.2.tar.gz:
	echo "obtain nlopt library"
	$(call DOWNLOAD,http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz)

rapidjson-1.0.2.tar.gz:
	echo "obtain rapidjson"
	$(call DOWNLOAD_RENAME,https://github.com/miloyip/rapidjson/archive/v1.0.2.tar.gz,$@)

sqlite-autoconf-3200100.tar.gz:
	echo "obtain sqlite"
	$(call DOWNLOAD,https://sqlite.org/2017/sqlite-autoconf-3200100.tar.gz)

zstd-1.3.1.tar.gz:
	echo "obtain zstd"
	$(call DOWNLOAD_RENAME,https://github.com/facebook/zstd/archive/v1.3.1.tar.gz,$@)

cnpy.zip:
	echo "obtain cnpy"
	$(call DOWNLOAD_RENAME,https://github.com/rogersce/cnpy/archive/master.zip,$@)

libdeflate-1.3.tar.gz:
	echo "obtain cnpy"
	$(call DOWNLOAD_RENAME,https://github.com/ebiggers/libdeflate/archive/v1.3.tar.gz,$@)

htslib-1.9.tar.gz:
	echo "obtain htslib"
	$(call DOWNLOAD_RENAME,https://github.com/samtools/htslib/archive/1.9.tar.gz,$@)


clean:
	-rm -rf $(LIB)
deepclean:
	-rm -rf $(LIB)
	-rm -rf pcre-8.36
	-rm -rf gsl-1.16
	-rm -rf eigen-3.3.7
	-rm -rf bzip2-1.0.6
	-rm -rf rapidjson-1.0.2
	-rm -rf nlopt-2.4.2
	-rm -rf zlib-1.2.8
	-rm -rf sqlite-autoconf-3200100
	-rm -rf zstd-1.3.1
	-rm -rf cnpy


##################################################
# Retired library
##################################################
tabix: tabix-0.2.6.tar.bz2 zlib
	-mv tabix-0.2.6.tar.bz2* tabix-0.2.6.tar.bz2
	tar jxf tabix-0.2.6.tar.bz2
	ln -s -f tabix-0.2.6 tabix
	sed -i.bak 's:^CFLAGS=:CFLAGS=-I../zlib-1.2.8:;s:^LIBPATH=:LIBPATH=-L../zlib-1.2.8:'  tabix/Makefile
	(cd tabix-0.2.6; make)

samtools: samtools-0.1.19.tar.bz2 zlib
	tar jxf $<
	ln -s -f samtools-0.1.19 samtools
	sed -i.bak 's:^CFLAGS=:CFLAGS=-I../zlib-1.2.8 -I../../zlib-1.2.8:;s:^LIBPATH=:LIBPATH=-L../zlib-1.2.8 -L../../zlib-1.2.8:' samtools/Makefile samtools/bcftools/Makefile
	(cd samtools; make lib-recur)

