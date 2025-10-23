#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tue Oct 21 06:21:45 PM CEST 2025
@author: sarkarjj
@hostname: hydra2

DESCRIPTION
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import os, glob

project_path= os.path.abspath('..')
savepath= os.path.join(project_path,'reports/savefig.png')
flat_filename= glob.glob(os.path.join(project_path, 'products/flat*NB02.fits'))[0]
files= sorted(glob.glob('/home/sarkarjj/Downloads/data/*14*/*NB02.fits'))
with fits.open(files[0]) as ref_hdu:
    ref_img= ref_hdu[0].data
with fits.open(files[-1]) as test_hdu:
    test_img= test_hdu[0].data
with fits.open(flat_filename) as flat_hdu:
    flat= flat_hdu[0].data

ref_corr= ref_img/flat
test_corr= test_img/flat

VMN, VMX= 0, 22e3
row, col, size=2400, 1600, 200 
plt.figure()
plt.subplot(1,2,1)
plt.title('raw')
plt.imshow(test_img[row:row+size, col:col+size], vmin= VMN, vmax= VMX)
plt.subplot(1,2,2)
plt.title('corrected')
plt.imshow(test_corr[row:row+size, col:col+size], vmin= VMN, vmax= VMX)
plt.savefig(savepath, dpi=300)
plt.close()

plt.figure()
plt.subplot(1,2,1)
plt.title('raw')
plt.imshow(test_img, vmin= VMN, vmax= VMX)
plt.subplot(1,2,2)
plt.title('corrected')
plt.imshow(test_corr, vmin= VMN, vmax= VMX)
plt.show()
