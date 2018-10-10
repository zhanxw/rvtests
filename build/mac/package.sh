#!/bin/bash

set -xue

DEPLOY_FILE=dist/test_${PLATFORM}.tar.gz
#_@VERSION@_${BATTLESHIP_PLATFORM}.exe
#cp executable/main $DEPLOY_FILE
tar zvcf $DEPLOY_FILE executable/main

set +xue
