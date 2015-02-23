#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Interactive plot to help create a Matplotlib colormap
(of the quantative kind, with ~linearly increasing luminance)

Pierre Haessig — February 2015
"""

from __future__ import division, print_function, unicode_literals

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from skimage import color


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
    

def lab2rgb_list(lab):
    '''wrapper of `lab2rgb` for list of tuples and (N,3) arrays'''
    lab = np.asarray(lab, dtype=float)
    assert lab.ndim == 2
    return color.lab2rgb(lab[None,:,:])[0,:,:]


def plot_cmap(cmap, **ax_kwargs):
    x = np.linspace(0,1, 500)[None, :]
    rgba = cmap(x)    

    fig, ax = plt.subplots(1,1, figsize=(6,2), num='Colormap demo')

    ax.imshow(rgba, interpolation="nearest")
    ax.set_aspect('auto')
    # remove ticks and use a light gray frame
    ax.set_xticks([])
    ax.set_yticks([])
    plt.setp(ax.spines.values(), color=(0.8,)*3)
    ax.set(**ax_kwargs)
    return ax

# blue green yellow
pts_lab = [
    (25, 50, -70),
    (50, -50, 50),
    (75, 20, 75),
]    

def plot_lab_pts(pts_lab, **plotargs):
    
    fig = plt.figure('Colormap creator')
    ax = fig.add_subplot(111)
    ax.set(
        xlabel = 'a', ylabel = 'b',
        xlim = (a_min,a_max),
        ylim = (b_min,b_max),
        aspect='equal',
    )

    # compute the sRGB color
    pts_rgb = lab2rgb_list(pts_lab)
    
    # Plot each point as an individual Line object
    for i in range(len(pts_lab)):
        l, a, b = pts_lab[i]
        plt.plot(a, b, 'o', color=pts_rgb[i], **plotargs)
    
    return ax
    
    
