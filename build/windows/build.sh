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

export CC=${MXE_DIR}/usr/bin/${MXE_TARGET}-gcc
export CXX=${MXE_DIR}/usr/bin/${MXE_TARGET}-g++
export LD=${MXE_DIR}/usr/bin/${MXE_TARGET}-ld
export AR=${MXE_DIR}/usr/bin/${MXE_TARGET}-ar
export FC=${MXE_DIR}/usr/bin/${MXE_TARGET}-gfortran
export PKG_CONFIG=${MXE_DIR}/usr/bin/${MXE_TARGET}-pkg-config
export LDFLAGS=

if [ "$PLATFORM" = "windows32" ]; then
    make -f Makefile.win CROSS=i686-w64-mingw32.static WIN32=1
fi

if [ "$PLATFORM" = "windows64" ]; then
    make -f Makefile.win CROSS=x86_64-w64-mingw32.static WIN32=1
fi

set +xue
