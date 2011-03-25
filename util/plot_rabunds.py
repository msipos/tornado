#!/usr/bin/python2

import pylab, numpy, scipy, sys

have_maxpy = False
try:
    import maxpy
    have_maxpy = True
except ImportError, e:
    pass

for f in sys.argv[1:]:
    d = numpy.loadtxt(f)

    if have_maxpy:
        dx = numpy.array(range(1, len(d)+1))
        x, y, e = maxpy.smart_bin(dx, d, 40)
        pylab.semilogy(x, y, '-x', label=f)
    else:
        pylab.semilogy(d, label=f)
pylab.legend(loc='upper right')
pylab.show()
