#!/bin/bash

set -xue

sudo apt-get update

# echo "deb http://pkg.mxe.cc/repos/apt/debian wheezy main" \
#         | sudo tee /etc/apt/sources.list.d/mxeapt.list
# sudo apt-key adv --keyserver x-hkp://keys.gnupg.net \
#         --recv-keys D43A795B73B16ABE9643FE1AFD8FFF16DB45C6AB

## pkg.mxe.cc is broken, now use a tweak
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 86B72ED9
#sudo add-apt-repository 'deb [arch=amd64] http://mirror.mxe.cc/repos/apt trusty main'
# or https for xenial
sudo add-apt-repository 'deb [arch=amd64] https://mirror.mxe.cc/repos/apt xenial main'
sudo apt-get update
apt-cache search mxe-

sudo apt-get update

sudo apt-get --yes install upx-ucl

if [ "$PLATFORM" = "windows32" ]; then
    MXE_TARGET=i686-w64-mingw32.static
fi

if [ "$PLATFORM" = "windows64" ]; then
    MXE_TARGET=x86-64-w64-mingw32.static
fi

MXE2_TARGET=$(echo "$MXE_TARGET" | sed 's/_/-/g')
sudo apt-get --yes install \
     mxe-${MXE2_TARGET}-qt

# MXE doesn't have 64bit NSIS
sudo apt-get --yes install \
     mxe-i686-w64-mingw32.static-nsis

sudo apt-get --yes install upx-ucl

set +xue
