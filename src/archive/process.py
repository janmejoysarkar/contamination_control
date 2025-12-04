#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tue Oct 21 06:21:45 PM CEST 2025
@author: sarkarjj
@hostname: hydra2

DESCRIPTION
- Runs the contamination correction on an image of choice.
- Saves the FITS file.
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import os, glob
import numpy as np

def calibration(filepath):
    # Open test image
    with fits.open(filepath) as test_hdu:
        test_img= test_hdu[0].data
        test_img_header= test_hdu[0].header
    FILT_NAME= test_img_header['FTR_NAME']

    # Open corresponding flat file
    if os.path.exists(flat_filename):
        with fits.open(flat_filename) as flat_hdu:
            flat= flat_hdu[0].data
    else:
        print('Check flat filepath')

    # Apply correction
    test_corr= test_img/flat
    if SAVE:
        save_fits_path=os.path.join(savepath, test_img_header['F_NAME'])
        fits.writeto(save_fits_path, test_corr, header=test_img_header, overwrite=True)
    return(test_img, test_corr)

def plotting(test_img, test_corr):
    plt.figure(figsize=(10,4))
    plt.subplot(1,2,1)
    plt.title('raw')
    plt.imshow(test_img[row-s:row+s, col-s:col+s], vmin= VMN)
    plt.colorbar()
    plt.subplot(1,2,2)
    plt.title('corrected')
    plt.imshow(test_corr[row-s:row+s, col-s:col+s], vmin= VMN)
    plt.colorbar()
    plt.savefig(preview_savepath, dpi=300)
    plt.close()

if __name__=='__main__':
    SAVE= True # Save corrected fits?
    SAVE_PRV= True # Save preview?
    # Paths
    project_path= os.path.abspath('..')
    filt_names=['NB01', 'NB02', 'NB03', 'NB04', 'NB05', 'NB06', 'NB07', 'NB08', 'BB01', 'BB02', 'BB03']
    
    # Run
    for filt_name in filt_names:
        savepath= os.path.join(project_path, 'products') 
        files= sorted(glob.glob(os.path.join(project_path, f'data/raw/*{filt_name}*'))) # Filepath for full disk images
        raw_filename= files[2]
        flat_filename= os.path.join(project_path, f'data/processed/flat_{os.path.basename(files[0])}')
        print(filt_name)
        print(os.path.basename(raw_filename), os.path.basename(flat_filename))
        test_image, corrected_image= calibration(raw_filename)
        # Plot if necessary
        if SAVE_PRV:
            VMN= 0
            row, col, s=1567, 1705, 200
            preview_savepath= os.path.join(project_path,f'reports/images/savefig_{os.path.basename(raw_filename)[:-5]}.png')
            plotting(test_image, corrected_image)

