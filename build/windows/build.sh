#!/bin/bash

set -xue

MXE_DIR=/usr/lib/mxe
export PATH=${MXE_DIR}/usr/bin:$PATH

if [ "$PLATFORM" = "windows32" ]; then
    MXE_TARGET=i686-w64-mingw32.static
fi

if [ "$PLATFORM" = "windows64" ]; then
    MXE_TARGET=x86_64-w64-mingw32.static
fi

export CC=${MXE_DIR}/usr/bin/i686-w64-mingw32.static-gcc
export CXX=${MXE_DIR}/usr/bin/i686-w64-mingw32.static-g++
export LD=${MXE_DIR}/usr/bin/i686-w64-mingw32.static-ld
export AR=${MXE_DIR}/usr/bin/i686-w64-mingw32.static-ar
export FC=${MXE_DIR}/usr/bin/i686-w64-mingw32.static-gfortran
export PKG_CONFIG=${MXE_DIR}/usr/bin/i686-w64-mingw32.static-pkg-config
make -f Makefile.win CROSS=i686-w64-mingw32.static- WIN32=1

set +xue
