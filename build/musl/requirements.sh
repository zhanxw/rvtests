#!/bin/bash

set -xue

sudo apt-get update

sudo apt-get --yes install \
    g++ make cmake \
    gfortran \
    wget
    #libqt4-dev upx-ucl \

wget https://skarnet.org/toolchains/native/x86_64-linux-musl-8.2.0.tar.xz
pwd
tar xf x86_64-linux-musl-8.2.0.tar.xz

set +xue
