import numpy, pylab

data = numpy.loadtxt('benchmark-results')
ns = data[:,0]
clink = data[:,3]
mothur = data[:,6]

pylab.loglog(ns, clink, 'rx', label = 'c-linkage')
pylab.loglog(ns, (ns ** 2) / 3e5, 'r--', label = 'N^2')
pylab.loglog(ns, mothur, 'bo', label = 'mothur')
pylab.loglog(ns, (ns ** 3) / 1e8, 'b--', label = 'N^3')
pylab.xlabel('N')
pylab.ylabel('time (seconds)')
pylab.axis([800, 10000, 1,  15000])
pylab.legend(loc = 'upper left')
pylab.show()
