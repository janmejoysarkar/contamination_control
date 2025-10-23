#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2025-07-21 11:31:32
@author: janmejoyarch
@hostname: suitpoc1

DESCRIPTION
"""

import numpy as np
import sunpy.map
from sunkit_image.coalignment import calculate_match_template_shift
from sunkit_image.coalignment import mapsequence_coalign_by_match_template
import os, glob
import matplotlib.pyplot as plt
import astropy.units as u


PLOT=True
project_path= os.path.abspath('..')
files= sorted(glob.glob(os.path.join(project_path, 'data/raw/set2/*')))
maps= sunpy.map.Map(files[1:], sequence=True)
ref_map= sunpy.map.Map(files[0])
#template
x1,y1,x2,y2= 293, 207, 483, 373
bottom_left = (min(x1, x2), min(y1, y2)) * u.pix
top_right = (max(x1, x2), max(y1, y2)) * u.pix
submap = ref_map.submap(bottom_left=bottom_left, top_right=top_right)
maps_aligned= mapsequence_coalign_by_match_template(maps, layer_index=0, clip=False, template=submap)
ani= maps_aligned.plot()
plt.show()

ff= np.ones(shape=(704,704))

for i in np.arange(1,len(maps)):
    m=maps[i].data
    s= m/ff
    M= maps_aligned[i].data
    ff= M/s

if PLOT:
    plt.figure()
    plt.subplot(1,3,1)
    plt.imshow(m, origin='lower')
    plt.subplot(1,3,2)
    plt.imshow(M, origin='lower')
    plt.subplot(1,3,3)
    plt.imshow(ff, origin='lower', vmin=0, vmax=2)
    plt.show()
