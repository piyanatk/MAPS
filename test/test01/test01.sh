#!/bin/sh

maps2uvfits test01.vis test01.uvfits -26.4311 117.357 400.0 $SIM/array/mwa_32.txt
rm -rf test01.uv
fits in='test01.uvfits' out='test01.uv' op='uvin'
rm -rf test01.map test01.beam test01.clean test01.restor
# Imaging both IFs at once in Miriad seems to change the channel bandwidth to be (max-min)/nchan
# which obviously only works for a continuous band, and not if there are separate IFs, so need to
# image separately, which can be done by using either of the below options in invert:
# 1)  If you want to select the frequencies based on spectral window (IF), select=window\(1\)
# 2)  If you want to select the frequencies based on channel number, line=channel,4,5
invert vis=test01.uv beam=test01.beam map=test01.map imsize=256,256
clean map=test01.map beam=test01.beam out=test01.clean niters=1000
restor model=test01.clean beam=test01.beam map=test01.map out=test01.restor
rm -rf test01.fits
fits op='xyout' in='test01.restor' out='test01.fits'
