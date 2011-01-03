#!/bin/bash

time mothur "#read.dist(column=$1.dist,name=$1.names); cluster(cutoff=0.1);" > /dev/null
