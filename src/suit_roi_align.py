#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mon Nov  3 05:34:31 PM CET 2025
@author: janmejoy
@hostname: machine

DESCRIPTION
"""

import matplotlib.pyplot as plt 
import glob 
from sunpy.coordinates import Helioprojective, SphericalScreen, propagate_with_solar_surface
from sunpy.map import Map, MapSequence
from astropy.coordinates import SkyCoord
import astropy.units as u
from sunkit_image.coalignment import apply_shifts, calculate_match_template_shift
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign
import os
from matplotlib.widgets import RectangleSelector
import numpy as np
from astropy.io import fits

def select_roi_with_mouse(sunpy_map, cmap=None, norm=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=sunpy_map)
    ax.set_title("Select ROI (click and drag) then close the window")
    sunpy_map.plot(axes=ax)
    coords = []
    def onselect(eclick, erelease):
        coords.append((eclick.xdata, eclick.ydata, erelease.xdata, erelease.ydata))

    toggle_selector = RectangleSelector(ax, onselect, useblit=True,
                      button=[1], minspanx=5, minspany=5, spancoords='pixels',
                      interactive=True)
    plt.show()

    if not coords:
        raise RuntimeError("ROI selection cancelled or failed.")

    x1, y1, x2, y2 = coords[0]
    
    bottom_left = (min(x1, x2), min(y1, y2)) * u.pix
    top_right = (max(x1, x2), max(y1, y2)) * u.pix

    submap = sunpy_map.submap(bottom_left=bottom_left, top_right=top_right)
    return submap

def align_maps(seq):
    '''
    Align a map sequence to the reference frame
    
    RETURNS:
    Aligned sunpy map sequence
    '''
    # Ref submap is taken from the first frame of the sequence
    ref_submap = select_roi_with_mouse(seq[0]) 
    nt = len(seq)
    xshift_keep = np.zeros(nt) * u.pix
    yshift_keep = np.zeros_like(xshift_keep)
    shifts = calculate_match_template_shift(seq, template=ref_submap,)
    xshift_arcseconds = shifts["x"]
    yshift_arcseconds = shifts["y"]
    for i, m in enumerate(seq):
        xshift_keep[i] = xshift_arcseconds[i] / m.scale[0]
        yshift_keep[i] = yshift_arcseconds[i] / m.scale[1]
    seq = apply_shifts(seq, -yshift_keep, -xshift_keep, clip=True)
    final_seq_2 = []
    for i,j in enumerate(seq):
        date = j.date.strftime('%H:%M:%S')
        dhobt_dt = j.meta['dhobt_dt']
        grt_dt = j.meta['grt_dt']
        mfgdate = j.meta['mfgdate']
        t_obs = j.meta['t_obs']
        date_obs = j.meta['date-obs']
        obs_strt = j.meta['obs_strt']
        obs_stp = j.meta['obs_stp']
        crtime = j.meta['crtime']
        exptime = j.meta['cmd_expt']
        meas_exptime = j.meta['meas_exp']
        p = Map(j.data, seq[0].meta)
        p.meta['dhobt_dt'] = dhobt_dt
        p.meta['grt_dt'] = grt_dt
        p.meta['mfgdate'] = mfgdate
        p.meta['t_obs'] = t_obs
        p.meta['date-obs'] = date_obs
        p.meta['obs_strt'] = obs_strt
        p.meta['obs_stp'] = obs_stp
        p.meta['crtime'] = crtime
        p.meta['cmd_expt'] = exptime
        p.meta['meas_exp'] = meas_exptime
        final_seq_2.append(p)
    final_seq_2 = Map(final_seq_2, sequence=True)
    return (final_seq_2)

def generate_flat(final_seq_2, SAVE=False):
    '''
    Generate a flat field of contaminant spots.
    RETURNS:
    Flat field file. Give the option to save or not.
    '''
    ref_map=final_seq_2[0]
    aligned_maps= [m.data for m in final_seq_2]

    if (MODE=='median'):
        aligned_map_arr= np.stack(aligned_maps, axis=0)
        combined_image= np.median(aligned_map_arr, axis=0)
    elif (MODE=='max'):
        aligned_map_arr= np.stack(aligned_maps, axis=0)
        combined_image= np.max(aligned_map_arr, axis=0)
    else:
        print("Specify image combination mode")
    flat= ref_map.data/combined_image
    if SAVE:
        fits.writeto(savepath, flat, overwrite= True)
        print('Flat file saved at', savepath)
    return flat

if __name__=='__main__':
    MODE='median'
    project_path = os.path.abspath("..")
    savepath= os.path.join(project_path, "data/interim/flat_frame.fits")
    f_seq = sorted(glob.glob(os.path.join(project_path,'data/raw/*')))
    seq = Map(f_seq, sequence=True)
    aligned_sequence= align_maps(seq)
    flat_frame= generate_flat(aligned_sequence, SAVE=True)
    
    ref_map= aligned_sequence[0]
    corrected_map= Map(aligned_sequence[0].data/flat_frame, ref_map.meta)
    corrected_map_ls= [Map(m.data/flat_frame, m.meta) for m in aligned_sequence]
    corrected_map_seq= MapSequence(corrected_map_ls)
    
    '''
    VMN= 0
    VMX= np.max(ref_map.data)
    fig, ax= plt.subplots(1,3, sharex=True, sharey=True)
    im0= ax[0].imshow(ref_map.data, origin='lower', vmin=VMN)
    ax[0].set_title('Raw')
    plt.colorbar(im0, ax=ax[0])
    
    im1= ax[1].imshow(flat, origin='lower', vmin=VMN, vmax=1.2)
    ax[1].set_title('Calibration frame')
    plt.colorbar(im1, ax=ax[1])
    
    im2= ax[2].imshow(corrected_img, origin='lower', vmin=VMN, vmax=VMX)
    ax[2].set_title('Corrected map')
    plt.colorbar(im2, ax=ax[2])
    plt.show()

    #save_fits
    save_path= os.path.join(project_path, 'data/processed/savefits.fits')
    with fits.open(f_seq[0]) as hdu:
        head= hdu[0].header
    fits.writeto(save_path, corrected_img, header= head, overwrite= True)
    print('File saved at', savepath)
    '''
