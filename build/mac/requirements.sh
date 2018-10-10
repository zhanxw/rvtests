#!/bin/bash

set -xue

brew update && brew install gcc; brew reinstall gcc --without-multilib
# brew install cmake

set +xue
