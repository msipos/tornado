#!/bin/bash

# Debugging flags
#CFLAGS="-Wall -g"

# Release flags
CFLAGS="-Wall -DNDEBUG -O3 -funroll-loops"

# LIBS
LIBS="-lm"

g++ $CFLAGS $LIBS -o tree-distance tree-distance.cpp

