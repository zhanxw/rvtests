#!/bin/bash

set -xue

sudo apt-get update

sudo apt-get --yes install \
    g++ make cmake \
    #libqt4-dev upx-ucl \

set +xue
