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

PLOT= True
SAVE= True
project_path= os.path.abspath('..')
savepath= os.path.join(project_path,'reports/savefig.png')
files= sorted(glob.glob('/home/sarkarjj/data/raw/*.fits'))
with fits.open(files[0]) as ref_hdu:
    ref_img= ref_hdu[0].data
    ref_header= ref_hdu[0].header
with fits.open(files[-1]) as test_hdu:
    test_img= test_hdu[0].data
FILT_NAME= ref_header['FTR_NAME']
flat_filename= glob.glob(os.path.join(project_path, f'data/processed/flat*{FILT_NAME}.fits'))[0]
with fits.open(flat_filename) as flat_hdu:
    flat= flat_hdu[0].data

ref_corr= ref_img/flat
test_corr= test_img/flat

VMN= 0
row, col, size=2400, 1600, 200 

if SAVE:
    plt.figure()
    plt.subplot(1,2,1)
    plt.title('raw')
    plt.imshow(test_img[row:row+size, col:col+size], vmin= VMN)
    plt.subplot(1,2,2)
    plt.title('corrected')
    plt.imshow(test_corr[row:row+size, col:col+size], vmin= VMN)
    plt.savefig(savepath, dpi=300)
    plt.close()

if PLOT:
    plt.figure()
    plt.subplot(1,2,1)
    plt.title('raw')
    plt.imshow(test_img, vmin= VMN)
    plt.subplot(1,2,2)
    plt.title('corrected')
    plt.imshow(test_corr, vmin= VMN)
    plt.show()
