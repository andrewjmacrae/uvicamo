# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import meep as mp
import numpy as np
import matplotlib.pyplot as plt

#This is a simulation of a field propagating through a waveguide bend
cell = mp.Vector3(16,16,0)
geometry = [mp.Block(mp.Vector3(12,1,mp.inf),
                     center=mp.Vector3(-2.5,-3.5),
                     material=mp.Medium(epsilon=12)),
            mp.Block(mp.Vector3(1,12,mp.inf),
                     center=mp.Vector3(3.5,2),
                     material=mp.Medium(epsilon=12))]
pml_layers = [mp.PML(1.0)]
resolution = 10

# we have two blocks to produce the bent waveguide

sources = [mp.Source(mp.ContinuousSource(wavelength=2*(11**0.5), width=20),
                     component=mp.Ez,
                     center=mp.Vector3(-7,-3.5),
                     size=mp.Vector3(0,1))]
#uses width = 20 so that it does not just turn on source at t = 0 since it 
#excites other frequencies so source is turned up slowly with hyperbolic tangent 
#function done automatically through meep

# the size object in sources no longer makes point source, it is line source
#with same width as waveguide

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution)

sim.run(until=200)

#Instead of running output_efield_z only at the end of the simulation
#we run it at every 0.6 time units via mp.at_every(0.6, mp.output_efield_z)
#%%

eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.axis('off')
plt.show()
#%%
ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.show()