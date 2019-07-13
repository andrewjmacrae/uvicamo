#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 16:03:27 2019

@author: noahgorgichuk
"""

import meep as mp
import numpy as np 
import matplotlib.pyplot as plt
from meep.materials import Cu

wid = 1
sq1 = 10
sq2 = 20
cell = mp.Vector3(2 * sq2, 2 * sq2, 0)
gap = 2

cell = cell = mp.Vector3(1.5 * sq2,1.5 * sq2,0)
#first square

#vertial
sq1_1 = mp.Block(mp.Vector3(sq1,wid,mp.inf), center=mp.Vector3(0,-sq1/2), material= mp.Medium(epsilon = 10))
sq1_2 = mp.Block(mp.Vector3(sq1,wid,mp.inf), center=mp.Vector3(0,sq1/2), material= mp.Medium(epsilon = 10))
#horizontal
sq1_3 = mp.Block(mp.Vector3(wid,sq1+wid,mp.inf), center=mp.Vector3(sq1/2,0), material= mp.Medium(epsilon = 10))
#sides from the split on first square

#length of split side
split_1 = abs(sq1/2 - gap/2)

sq1_4 = mp.Block(mp.Vector3(wid,(sq1/2) - (gap / 2),mp.inf), center=mp.Vector3(-sq1/2,(sq1 - split_1 + wid)/2), material= mp.Medium(epsilon = 10))
sq1_5 = mp.Block(mp.Vector3(wid,(sq1/2) - (gap / 2),mp.inf), center=mp.Vector3(-sq1/2,(-sq1 + split_1 - wid)/2), material= mp.Medium(epsilon = 10))

#second square

#vertial
sq2_1 = mp.Block(mp.Vector3(sq2,wid,mp.inf), center=mp.Vector3(0,-sq2/2), material= mp.Medium(epsilon = 10))
sq2_2 = mp.Block(mp.Vector3(sq2,wid,mp.inf), center=mp.Vector3(0,sq2/2), material= mp.Medium(epsilon = 10))
#horizontal
sq2_3 = mp.Block(mp.Vector3(wid,sq2+wid,mp.inf), center=mp.Vector3(-sq2/2,0), material= mp.Medium(epsilon = 10))
#sides from the split on first square
split_2 = abs(sq2/2 - gap/2)
sq2_4 = mp.Block(mp.Vector3(wid,(sq2/2) - (gap / 2),mp.inf), center=mp.Vector3(sq2/2,(-sq2 + split_2 - wid)/2), material= mp.Medium(epsilon = 10))
sq2_5 = mp.Block(mp.Vector3(wid,(sq2/2) - (gap / 2),mp.inf), center=mp.Vector3(sq2/2,(sq2 - split_2 + wid)/2), material= mp.Medium(epsilon = 10))


pml_layers = [mp.PML(1.0)]

resolution = 10

geometry = [sq1_1,sq1_2,sq1_3,sq1_4,sq1_5,sq2_1,sq2_2,sq2_3,sq2_4,sq2_5]

#10.3GHz frequency has a wavelength of 29106 um
sources = [mp.Source(mp.ContinuousSource(wavelength=10, width = 20),
                     component=mp.Ez,
                     center=mp.Vector3(sq1,0),
                     size=mp.Vector3(0,1))]


sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution)
sim.run(until = 200)

#%%
eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
#get_array obtains a slice of data which outputs to NumPy array and displays the results

plt.figure()
plt.imshow(eps_data, interpolation='spline36', cmap='binary')
plt.axis('on')
plt.grid('on')
plt.show()
#%%
ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.show()