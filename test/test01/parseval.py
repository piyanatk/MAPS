#!/usr/bin/python
#
# Check the Parseval's identity for MAPS_im2uv program
# from the MAPS package.
# For the source brightness image, B_{i,j}, and the output
# UV-visibility image, V_{k,l}, the Parseval's identity states that
#
#    \sum_{i,j}^{N^2} |B_{i,j}|^2=\sum_{k,l}^{N^2} |V_{k,l}|^2,
#
# i.e. "the total power in the brightness domain is equal to the total
# power in its FFT image (here 2-dimensional).
#
# This is to help find correct normalization factor to use in MAPS_im2uv.
#
# The file 'brimout.txt' is written by LOsim.
# The file 'test01_Visibility.dat' is written by MAPS_im2uv.
#

from pylab import *

b = loadtxt('brimout.txt') # Brightness

v = fromfile('test01_Visibility.dat', double) # Complex visibilities
v = v[1:] # First two int32 numbers, the size, occupy 64 bits of 1 double 
v = v.reshape((1024,2048))

rev = v[:,0::2] # Real parts: every other number starting with the first
imv = v[:,1::2] # Imag parts: every other number starting with the second

v2 = rev**2 + imv**2  # \sum |V|^2 as complex numbers

print 'Sum(B^2) = ', sum(b**2), \
      ', Sum(V^2) = ', sum(v2)
#print 'Sum(B)/1024**2 = ', sum(b**2)/1024**2

#
# Plot brightness and visibility images
#
Vcp = rev + 1j*imv # Make complex visibility in V
Vsh = fftshift(Vcp)
V2 = real(Vsh*conj(Vsh))
V = sqrt(V2)

figure(figsize=[12,4])
subplot(131)
imshow(b, cmap=cm.hot); xlabel('Brightness image')
subplot(132)
imshow(sqrt(v2), cmap=cm.hot); xlabel('Visibility image')
subplot(133)
imshow(sqrt(V2), cmap=cm.hot); xlabel('Visibility fft-shifted')

show()
