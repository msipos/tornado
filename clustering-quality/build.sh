#!/bin/bash

CFLAGS="`pkg-config --cflags glib-2.0` -std=gnu99"
LIBS=`pkg-config --libs glib-2.0`

gcc $LIBS $CFLAGS cluster-quality.c cfasta.c -o cluster-quality
