#!/bin/bash

set -xue

brew update && brew install gcc && brew link gcc

set +xue
