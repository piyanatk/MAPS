#!/usr/bin/python
#
# Normalyzers for N-dimensional FFT 
#
# The conclusion is the normalyzer must be 1/N^d,
# where d is number of dimensions, to conform with
# the Parseval's identity
#

from pylab import *

N = 100

x = rand(N)
y = fft(x)
y2 = abs(y)**2
print 'N = 1: sum(x^2) = ', sum(x**2), ', sum(y^2) = ', sum(y2), \
      ', sum(y^2)/N = ', real(sum(y2))/float(N)  

x = rand(N,N)
y = fft2(x)
y2 = abs(y)**2
print 'N = 2: sum(x^2) = ', sum(x**2), ', sum(y^2) = ', sum(y2), \
      ', sum(y^2)/N^2 = ', real(sum(y2))/float(N**2)  

x = rand(N,N,N)
y = fftn(x)
y2 = abs(y)**2
print 'N = 3: sum(x^2) = ', sum(x**2), ', sum(y^2) = ', sum(y2), \
      ', sum(y^2)/N^3 = ', real(sum(y2))/float(N**3)  

N = 32

x = rand(N,N,N,N)
y = fftn(x)
y2 = abs(y)**2
print 'N = 4: sum(x^2) = ', sum(x**2), ', sum(y^2) = ', sum(y2), \
      ', sum(y^2)/N^4 = ', real(sum(y2))/float(N**4)  

N = 16

x = rand(N,N,N,N,N)
y = fftn(x)
y2 = abs(y)**2
print 'N = 5: sum(x^2) = ', sum(x**2), ', sum(y^2) = ', sum(y2), \
      ', sum(y^2)/N^5 = ', real(sum(y2))/float(N**5)  

x = rand(N,N,N,N,N,N)
y = fftn(x)
y2 = abs(y)**2
print 'N = 6: sum(x^2) = ', sum(x**2), ', sum(y^2) = ', sum(y2), \
      ', sum(y^2)/N^6 = ', real(sum(y2))/float(N**6)  

