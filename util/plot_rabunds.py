#!/usr/bin/python2

import pylab, numpy, scipy, sys

for f in sys.argv[1:]:
    d = numpy.loadtxt(f)
    pylab.semilogy(d, label=f)
pylab.legend(loc='upper right')
pylab.show()
