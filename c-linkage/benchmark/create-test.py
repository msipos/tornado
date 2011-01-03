import sys, random

name = sys.argv[1]
num = int(sys.argv[2])

maxdist = 0.2

nf = open('%s.names' % name, 'w')
for i in xrange(1, num+1):
  nf.write("%d %d\n" % (i, i))
nf.close()

df = open('%s.dist' % name, 'w')
for i in xrange(1, num+1):
  for j in xrange(i+1, num+1):
    df.write("%d %d %f\n" % (i, j, random.uniform(0.0, maxdist)))
df.close()
