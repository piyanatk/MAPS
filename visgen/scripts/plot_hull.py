#!/usr/bin/python
#
# Plot the convex hull of N points on a plane
#
from pylab import *

def plot_border(bd, r):
    nbd = len(bd)
    N = r.shape[0]
    for i in xrange(nbd):
        i1 = mod((i+1), nbd)   # If i+1 < nbd, i1 = i+1, if i == nbd, i1 = 0
        plot((r[bd[i],0], r[bd[i1],0]), (r[bd[i],1], r[bd[i1],1]), 'b')
        plot(r[bd[i],0], r[bd[i],1], 'wo')
        plot(r[bd[i],0], r[bd[i],1], 'g.')
        axis('equal'); grid(1)

r = loadtxt('points.txt')
bd = loadtxt('bd_idx.txt')

plot(r[:,0], r[:,1], 'g.'); axis('equal'); grid(1)
plot_border(bd, r)

show()
