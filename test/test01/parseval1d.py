#!/usr/bin/python
from pylab import *

x = rand(1024) + sin(pi*rand(1024)/510.) + \
    5*sin(6*pi*linspace(0,1,1024)+rand(1024)) + \
    4*sin(4.5*pi*linspace(0,1,1024)+rand(1024))

y = fftshift(fft(x))
y2 = y*conj(y)


print 'sum(x^2) = ', sum(x**2), ', sum(y^2) = ', real(sum(y2)), \
      'sum(y^2)/N = ', real(sum(y2))/1024.  

f = linspace(-512, 512, 1025)
y = hstack((y,y[0])) # y[0] is Nyquist freq for even N

figure(figsize=[9,4])
subplot(121); plot(x, 'b'); grid(1); xlabel('x(t)')
subplot(122); plot(f, abs(y), 'r'); grid(1); xlabel('FFT(x(t))')

show()

