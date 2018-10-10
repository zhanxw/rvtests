#!/bin/bash

set -xue

DEPLOY_FILE=dist/test_${PLATFORM}.tar.gz
#_@VERSION@_${BATTLESHIP_PLATFORM}.exe
#cp executable/main $DEPLOY_FILE
rm -f executable/*.d
tar zvcf $DEPLOY_FILE executable example README.md

set +xue
