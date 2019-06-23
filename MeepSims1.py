# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import meep as mp
#field pattern viewed with a localized Continuous Wave source in a waveguide
#source at one end and watch the fields propagate 
#down in x direction

#creates cell length of 16 um (x-direction), 8 um (y-direction) and 0 (z-direction)
#vector3 object stores the size of cell
cell = mp.Vector3(16,8,0)

#Geometric object stored in the geometry object:This class, and its descendants, 
#are used to specify the solid geometric objects that form the dielectric structure
# being simulated.


geometry = [mp.Block(mp.Vector3(mp.inf,1,mp.inf),
                     center=mp.Vector3(),
                     material=mp.Medium(epsilon=12))]

#block is the shape of material given by Vector 3
#The waveguide is specified by a Block (parallelepiped) of size ∞×1×∞, with ε=12
#centered at (0,0) which is the center of the cell. By default, any place where 
#there are no objects there is air (ε=1)

sources = [mp.Source(mp.ContinuousSource(frequency=0.15),
                     component=mp.Ez,
                     center=mp.Vector3(-7,0))]
#We gave the source a frequency of 0.15, and specified a ContinuousSource which is 
#just a fixed-frequency sinusoid exp(−iωt) that by default is turned on at t=0
# frequency is specified in units of 2πc, which is equivalent to the inverse of the 
#vacuum wavelength. Thus, 0.15 corresponds to a vacuum wavelength of about 1/0.15=6.67 μm, 
#or a wavelength of about 2 μm in the ε=12 material 

pml_layers = [mp.PML(1.0)]
#boundry conditions that absorb around the cell
#absorbing boundries handled by perfectly matched layers (PML)
#this adds absorbing layer thickness of 1 um around all sides of cell

resolution = 10
#this gives resolution of 10 pixels per um

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution)
#this is the simulation based on all previously defined objects

sim.run(until = 200)
#%%

import numpy as np
import matplotlib.pyplot as plt

eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
#get_array obtains a slice of data which outputs to NumPy array and displays the results

plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.axis('off')
plt.show()
#%%
#Next, we create an image of the scalar electric field Ez by 
#overlaying the dielectric function

ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.show()
#We see that the the source has excited the waveguide mode but has also excited 
#radiating fields propagating away from the waveguide. At the boundaries, the field 
#quickly goes to zero due to the PML.

#dark red is negative, white is zero, and dark blue is positive