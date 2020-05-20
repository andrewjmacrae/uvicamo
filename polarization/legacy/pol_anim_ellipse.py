# Animation of polarization ellipse for various angles due to imperfect WP

import numpy as np
from matplotlib import animation, pyplot as plt
pi = np.pi


# Definately should make this a module/class
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

# retardance of waveplate (ie lambda/ret)
ret = 4.523
# Initial Jones vector
j0 = (np.array([1,0]))+0j
j0/= np.sqrt(np.dot(j0,np.conjugate(j0))) # Normalize it!

#Define handles to plot for animation
fig, ax = plt.subplots()

# Set up plot parameters
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_title('$\lambda$/'+str(ret)+' Waveplate')
ax.set_aspect(1.0)
line, = ax.plot([],[])
txt = ax.text(-.9,.8,'')

def init_anim():
    return line,txt,

def run_anim(tht):
    global j0,ret
    wt = np.linspace(0,2*pi,180)
    j = np.dot(wp_r(2*pi/ret,tht),j0)
    Ex = j[0]*np.exp(1j*wt)
    Ey = j[1]*np.exp(1j*wt)
    line.set_data(np.real(Ex),np.real(Ey))
    txt.set_text('WP @'+str(round(180*tht/pi))+'$^\circ$')
    return line,txt,

anm = animation.FuncAnimation(fig,run_anim,interval = 60, frames = np.linspace(0,pi,100),init_func=init_anim,blit= True)
plt.show()