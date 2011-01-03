#!/bin/bash

# Debugging flags
#CFLAGS="-Wall -g -pg"
# Release flags
CFLAGS="-Wall -DNDEBUG -O3 -funroll-loops"

# LIBS
LIBS="-lm"

g++ $CFLAGS $LIBS -o c-linkage cluster.cpp load.cpp dist.cpp
