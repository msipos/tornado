#!/bin/bash

# Debugging flags
#CFLAGS="-Wall -pg -g"
#LIBS="-lgsl -lgslcblas -lm"

# Release flags
CFLAGS="-Wall -DNDEBUG -O3 -funroll-loops"
LIBS="-lm"

g++ $CFLAGS $LIBS -o psi-distance psi-distance.cpp
