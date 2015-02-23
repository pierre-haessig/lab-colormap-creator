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


def rgb_saturate(rgb, gray_out=0.2):
    '''inplace saturation of values outside [0,1].
    `rgb` should be a (..., 3) shaped array
    `gray_out` should be between 0 and 1
    
    '''
    rgb = np.atleast_2d(rgb)
    assert 0 <= gray_out <= 1
    # Detect pixel saturation, in *any* of the R,G,B channel
    sat = ((rgb<0) | (rgb>1)).sum(axis=-1).astype(bool)

    # Apply channel saturation
    rgb[rgb>1] = 1
    rgb[rgb<0] = 0

    # Gray out pixels in the saturated zones
    rgb[sat] = 0.5 + (rgb[sat]-0.5)*gray_out
    
    return rgb, sat

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
    ab_grid_rgb, _ = rgb_saturate(ab_grid_rgb)
    
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
    '''plot a matplotlib colormap, as an horizontal bar'''
    x = np.linspace(0,1, 500)[None, :]
    rgba = cmap(x)    

    fig, ax = plt.subplots(1,1, figsize=(6,1), num='Colormap demo')

    ax.imshow(rgba, interpolation="nearest")
    ax.set_aspect('auto')
    # remove ticks and use a light gray frame
    ax.set_xticks([])
    ax.set_yticks([])
    plt.setp(ax.spines.values(), color=(0.8,)*3)
    ax.set(**ax_kwargs)
    return ax


def interp_lab(pts_lab, n):
    '''interpolates `pts_lab` linearly in the Lab space,
    by adding `n` times as many points
    
    number of output points: (len(pts)-1)*n + 1)
    '''
    # unpack Lab channels
    l,a,b = np.asarray(pts_lab).T
    # Interpolate each channel linearly
    xp = np.linspace(0,1, len(l))
    x = np.linspace(0,1,  (len(l)-1)*n + 1)
    l_interp, a_interp, b_interp = [
        np.interp(x, xp, channel) for channel in [l,a,b]
    ]
    # repack the Lab channels
    pts_lab_interp = np.vstack((l_interp, a_interp, b_interp)).T
    return pts_lab_interp


def cmap_from_lab_pts(pts_lab, n_interp=50):
    'Creates a colormap from the interpolation (in Lab space) of `pts_lab`'
    pts_lab_interp = interp_lab(pts_lab, n_interp)
    # compute the sRGB color
    pts_rgb = lab2rgb_list(pts_lab_interp)
    pts_rgb, sat = rgb_saturate(pts_rgb, gray_out=1)
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Lab cmap', pts_rgb, len(pts_lab_interp))
    return cmap


def plot_lab_pts(pts_lab, gray_out=0.2, **plotargs):
    fig = plt.figure('Colormap creator')
    ax = fig.add_subplot(111)
    ax.set(
        xlabel = 'a', ylabel = 'b',
        xlim = (a_min,a_max),
        ylim = (b_min,b_max),
        aspect='equal',
    )
    
    # add a marker at neutral gray and use a gray background
    ax.patch.set_color((0.85, 0.85, 0.85))
    
    # compute the sRGB color
    pts_rgb = lab2rgb_list(pts_lab)
    pts_rgb, sat = rgb_saturate(pts_rgb, gray_out)
    
    # Plot each point as an individual Line object
    for i in range(len(pts_lab)):
        l, a, b = pts_lab[i]
        mec = 'r' if sat[i] else 'k'
        mew = 1. if sat[i] else 0.
        plt.plot(a, b, 'o', color=pts_rgb[i],
            markeredgecolor=mec,
            markeredgewidth=mew,
            **plotargs)
    
    return ax
    
def demo_cmap(pts_lab, title=''):
    plot_cmap(cmap_from_lab_pts(pts_lab), title=title)
    ax = plot_lab_pts(interp_lab(pts_lab, 5), markersize=10)
    ax = plot_lab_pts(pts_lab, markersize=20)
    ax.plot(0,0,'+', color='k')
    # TODO: add the cmap test
    
    
### Interactive Lab colormap editor ###

class LabEditor(object):
    '''
    Interactive Lab points editor to create a quantitative colormap,
    with linearly increasing lightness (L channel).
    
    Usage
    -----
    # Initial L-a-b points to start with.
    # Notice the choice of a linear increase in L 
    >>> pts_lab = [
         [ 10., 40,-52],
         [ 30., -5,-20],
         [ 50.,-42, 18],
         [ 70.,-68, 67],
         [ 90., -9, 87]]
    
    # Launch the editor (with an interactive backend like qt),
    # and do you edit by clicking on the points in the (a,b) space
    >>> %matplotlib qt # if needed
    >>> ed = LabEditor(pts_lab)
 
    
    # After moving the points, retrieve the modified L-a-b points at:
    >>> led.pts_lab
    array([[ 10.        ,  36.41129032, -49.87903226],
           [ 30.        ,   0.84677419, -29.55645161],
           [ 50.        , -31.33064516,   3.46774194],
           [ 70.        , -31.33064516,  52.58064516],
           [ 90.        ,  -8.46774194,  88.99193548]])
    '''
    def __init__(self, pts_lab):
        self.edit_ind = None
        pts_lab = np.asarray(pts_lab, dtype=float).copy()
        self.pts_lab = pts_lab
        ax = plot_lab_pts(pts_lab, markersize=15, picker=5)
        self.ax = ax
        self.fig = ax.figure
        #self.line = ax.plot(pts_lab[:,1], pts_lab[:,2])
        plt.show()
        
        # connect events
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.canvas.mpl_connect('motion_notify_event', self.onmotion)
        self.fig.canvas.mpl_connect('button_release_event', self.onrelease)

    
    def onpick(self, event):
        marker = event.artist
        ind = self.ax.lines.index(marker)
        self.edit_ind = ind
        
        marker.set_markeredgewidth(2)
        
        # Find the lightness and add the ab plane in the background
        light = self.pts_lab[ind,0]
        ab_grid_rgb = ab_plane(light, 300, 300)
        self.ab_img = self.ax.imshow(ab_grid_rgb, **ab_imparams)
        
        self.fig.canvas.draw()
    
    def onmotion(self, event):
        'motion of a point in the (a,b) space. Update position and color'
        ind = self.edit_ind
        if ind is None:
            return
        marker = self.ax.lines[ind]
        a_new, b_new = event.xdata, event.ydata
        
        marker.set_xdata([a_new])
        marker.set_ydata([b_new])
        
        self.pts_lab[ind,1] = a_new
        self.pts_lab[ind,2] = b_new
        
        # Compute the new color:
        rgb = lab2rgb_list(self.pts_lab[[ind]])[0]
        rgb, sat = rgb_saturate(rgb)
        rgb = rgb[0]
        sat = sat[0]
        marker.set_markeredgecolor('r' if sat else 'k')
        marker.set_color(rgb)
        
        # update
        self.fig.canvas.draw()
    
    def onrelease(self, event):
        ind = self.edit_ind
        self.edit_ind = None
        if ind is not None:
            marker = self.ax.lines[ind]
            marker.set_markeredgewidth(0.5)
            
            self.ab_img.remove()
            
            self.fig.canvas.draw()
            #print('new color points:')
            #print(self.pts_lab)

