#!/usr/bin/python

# Reformat .cluster and .list files to be more readable
# Take the first line of the file, remove the first 2 words
# (clustering distance and number of OTUs).
# Print the remainder of the data split by newlines (not spaces).

import sys

f = open(sys.argv[1])
line = f.readline().strip()

arr = line.split()
for e in arr[2:]:
    print(e)
