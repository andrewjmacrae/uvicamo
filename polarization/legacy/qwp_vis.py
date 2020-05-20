# Simulate a fast rotating polarizer with a slowly rotating waveplate of imperfect delay
# Plots the envelope of amplitudes

import numpy as np
from matplotlib import pyplot as plt
pi = np.pi

def rot(th):
    ''' returns a rotation matrix. th is the angle of rotation in radians'''
    return np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
def wp(phi):
    ''' Jones matrix for a wave plate with fast axis horizontal
        - phi is phase shift. 
        - Retardance is phi * lambda/2pi so a QWP has phi = pi/2 '''
    return np.array([[1,0],[0,np.exp(1j*phi)]])
def wp_r(phi,th):
    ''' Waveplate of phase delay phi, rotated by angle th wrt horizontal '''
    return np.dot(rot(-th),np.dot(wp(phi),rot(th)))
def pol():
    ''' Jones matrix for a polarizer '''
    return np.array([[1,0],[0,0]])
def pol_r(th):
    ''' Polarizer rotated by angle th '''
    return np.dot(rot(-th),np.dot(pol(),rot(th)))

# Inital Jones matrix - assume horizontal
j0 = np.array([1,0])

# Pick retardance here 2*pi/x for a lambda/x plate
ret = 4.3
retardance = 2*pi/ret 

# range of angles to scan the polarizer
pol_angle  = np.linspace(0,2*pi,100)
# range of angles over which to scan the waveplate
phi = np.linspace(0,pi,500)
vis = phi*0
cnt = 0

# Inelegant code ... for each WP angle, spin the polarizer and find the maximum and minimum amplitudes
for phi0 in phi:
    Imin = 1
    Imax = -1
    for theta in pol_angle:
        T = np.dot(np.dot(pol_r(theta),wp_r(retardance,phi0)),j0)
        I = abs(T[0])**2 + abs(T[1])**2
        Imax = max(I,Imax)
        Imin = min(I,Imin)
        
    vis[cnt] = Imax - Imin
    cnt += 1

# ... and plot it
plt.plot(phi/pi,vis)

plt.xlim([0,1])
plt.grid(True)
plt.ylim([0,1.1])
plt.xlabel('Waveplate angle [$\\theta/2\pi$]')
plt.ylabel('Fringe Visibility [$I_{max}-I_{min}$]')
plt.title('Visibility for a $\lambda/$'+str(ret)+' Plate')
plt.show()