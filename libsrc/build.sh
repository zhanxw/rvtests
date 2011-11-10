#!/bin/sh
g++ -g -O0 -c *.cpp
rm lib-goncalo.a
ar -cqs lib-goncalo.a *.o
