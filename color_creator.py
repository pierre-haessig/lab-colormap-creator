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


na, nb =(20, 21)
L = 55.
a_min,a_max = (-0.7,1.)
b_min, b_max = (-0.9,0.8)

lab_grid = np.zeros((nb,na,3))
lab_grid[:,:,0] = L
b_grid, a_grid = np.ogrid[b_min:b_max:1j*nb,
                          a_min:a_max:1j*na]
lab_grid[:,:,1] = a_grid
lab_grid[:,:,2] = b_grid

rgb_grid = color.lab2rgb(lab_grid)

plt.show()
