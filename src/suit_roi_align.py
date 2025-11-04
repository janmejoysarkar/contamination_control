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

if __name__=='__main__':
    project_path = os.path.abspath("..")
    f_bb3 = sorted(glob.glob(os.path.join(project_path,'data/raw/*')))[5:]
    bb3 = Map(f_bb3, sequence=True)

    ref_submap = select_roi_with_mouse(bb3[0])#, norm=ImageNormalize(stretch=LinearStretch()))
    nt = len(bb3)

    xshift_keep = np.zeros(nt) * u.pix
    yshift_keep = np.zeros_like(xshift_keep)

    shifts = calculate_match_template_shift(bb3, template=ref_submap,)
    xshift_arcseconds = shifts["x"]
    yshift_arcseconds = shifts["y"]

    for i, m in enumerate(bb3):
        xshift_keep[i] = xshift_arcseconds[i] / m.scale[0]
        yshift_keep[i] = yshift_arcseconds[i] / m.scale[1]

    bb3 = apply_shifts(bb3, -yshift_keep, -xshift_keep, clip=True)
    final_seq_2 = []
    for i,j in enumerate(bb3):
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
        p = Map(j.data, bb3[0].meta)
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
    ref_map=final_seq_2[0]
    aligned_map_arr= np.stack([m.data for m in final_seq_2[5:]], axis=0)
    med= np.median(aligned_map_arr, axis=0)
    flat= ref_map.data/med

    plt.figure()
    plt.subplot(1,3,1)
    plt.imshow(ref_map.data, origin='lower', vmin=0)
    plt.title('Raw')
    plt.colorbar()
    plt.subplot(1,3,2)
    plt.imshow(flat, origin='lower', vmin=0, vmax=1.2)
    plt.title('Calibration frame')
    plt.colorbar()
    plt.subplot(1,3,3)
    plt.imshow(ref_map.data/flat, origin='lower', vmin=0, vmax=1.5e4)
    plt.title('Corrected map')
    plt.colorbar()
    plt.show()
