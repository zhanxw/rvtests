#!/bin/bash

set -xue
MXE_DIR=$HOME/build/zhanxw/rvtests/x86_64-linux-musl-8.2.0
export PATH=${MXE_DIR}/bin:$PATH
MXE_TARGET=x86_64-linux-musl

export CC=${MXE_DIR}/bin/${MXE_TARGET}-gcc
export CXX=${MXE_DIR}/bin/${MXE_TARGET}-g++
export LD=${MXE_DIR}/bin/${MXE_TARGET}-ld
export AR=${MXE_DIR}/bin/${MXE_TARGET}-ar
export FC=${MXE_DIR}/bin/${MXE_TARGET}-gfortran
export PKG_CONFIG=${MXE_DIR}/bin/${MXE_TARGET}-pkg-config
export LDFLAGS=

MXE_DIR=$HOME/x86_64-linux-musl-8.2.0
ls $HOME
ls $MXE_DIR
ls ${MXE_DIR}/bin
${CC}

make -f Makefile.win CROSS=${MXE_TARGET} MUSL=1

set +xue
