#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Interactive plot to help create a Matplotlib colormap
(of the quantative kind, with ~linearly increasing luminance)

Pierre Haessig — February 2015
"""

from __future__ import division, print_function, unicode_literals

import numpy as np
import matplotlib.pyplot as plt
from skimage import color


# create the initial plot

fig = plt.figure('Colormap creator')
ax = fig.add_subplot(111)
ax.set(xlabel='a', ylabel='b')


x0 = [0,1,2]
y0 = [0,1,1.5]
c = ['r', 'g', 'b']

lines = []

for i in range(len(x0)):
    lines.append(
    ax.plot(x0[i], y0[i], 'D', color=c[i])
    )


# Global extent of the ab plane
a_min,a_max = (-100.,100.)
b_min, b_max = (-110.,100.)

ab_imparams = dict(origin='lower', extent=[a_min, a_max, b_min, b_max],
                   interpolation='nearest', cmap='gray')


def ab_plane(L, na, nb):
    '''Generate a sRBG image of the ab plane at luminance `L`
    
    a increases along the axis 1, while
    b increases along axis 0.
    '''
    # Generate the ab grid
    ab_grid_lab = np.zeros((nb,na,3))
    ab_grid_lab[:,:,0] = L
    b_grid, a_grid = np.ogrid[b_min:b_max:1j*nb,
                              a_min:a_max:1j*na]
    ab_grid_lab[:,:,1] = a_grid
    ab_grid_lab[:,:,2] = b_grid

    # Convert Lab -> sRBG
    ab_grid_rgb = color.lab2rgb(ab_grid_lab)

    # Detect saturation in the RGB space:
    sat = ((ab_grid_rgb<0) | (ab_grid_rgb>1)).sum(axis=-1).astype(bool)
    #plt.imshow(sat, **ab_imparams)

    # Apply channel saturation
    ab_grid_rgb[ab_grid_rgb>1] = 1
    ab_grid_rgb[ab_grid_rgb<0] = 0

    # Gray out the saturated zones
    ab_grid_rgb[sat] = 0.5 + (ab_grid_rgb[sat]-0.5)*0.2
    
    return ab_grid_rgb


def plot_ab_plane(L):
    '''plot the output of `ab_plane` generation function'''
    ab_grid_rgb = ab_plane(L, 300, 300)
    plt.imshow(ab_grid_rgb, **ab_imparams)
    # add a small cross at the center (neutral gray)
    plt.plot(0,0,'+', color='w' if L<50 else 'k')
    # annotations:
    ax = plt.gca()
    ax.set(
        xlabel = 'a', ylabel = 'b',
        xlim = (a_min,a_max),
        ylim = (b_min,b_max),
        title = 'Lab space cut at L={:.0f}'.format(L)
        )
    plt.tight_layout()
    

plt.show()
