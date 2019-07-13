# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import meep as mp
from meep.materials import Cu


r = 1 #radius of inner ring
w = 1 #width of waveguide
dpml = 2 #thickness of PML 
pad = 4 #padding between waveguide and edge of PML
cell_xy = 4*(r+w+dpml+pad)

#first ring
#simulation causing problems when using copper (material = Cu) instead of whats below
c1 = mp.Cylinder(radius = r + w, material = mp.Medium(index = 3.4))
c2 = mp.Cylinder(radius = r)

#second ring
c3 = mp.Cylinder(radius = r + w + dpml + pad)
c4 = mp. Cylinder(radius = 1 + r + w +dpml + pad)

fcen = 0.118 #pulse center freq
df = 0.01 #pulse frequency width
src = mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Ez, mp.Vector3(r+0.1))

#In order to extract the frequencies and lifetimes (which may be 
#infinite in a lossless system) with FDTD, the basic strategy is simple. 
#You set up the structure with Bloch-periodic and/or absorbing boundaries,
# depending on whether it is a periodic or open system. Then you excite the mode(s) 
#with a short pulse (broad bandwidth) from a current placed directly inside 
#the cavity/waveguide/whatever.

#Once the source is turned off, fields boucing around the system can be analyzed

sim = mp.Simulation(cell_size=mp.Vector3(cell_xy, cell_xy),
                    geometry=[c1, c2, c3, c4],
                    sources=[src],
                    resolution=10,                    
                    boundary_layers=[mp.PML(dpml)])

sim.run(mp.at_beginning(mp.output_epsilon),
        mp.after_sources(mp.Harminv(mp.Ez, mp.Vector3(r+0.1), fcen, df)),
        until_after_sources=300)
# we run simulation until some additional period after sources finished
# mp.after_sources(mp.Harmiv(... identifies frequencies and decay rates of modes that 
#were excited after the source finished

#at the end harmiv prints list of frequencies it found
#first column is the real part of ωn, expressed in our usual 2πc 
#second collumn is imaginary part where negative means exponential decay
#third is "Q" which is lifetime (number of optical periods for energy decay by exp(-2*pi))
#fourth is absoute value of amplitude and fifth is complex amp

sim.run(mp.at_every(1/fcen/20, mp.output_efield_z), until=1/fcen)
#output 20 field snapshots over a whole period 1/fcen by appending the command