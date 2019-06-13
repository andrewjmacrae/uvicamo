import numpy as np
from matplotlib import pyplot as plt
pi = np.pi

def rot(th):
    return np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
def wp(phi):
    return np.array([[1,0],[0,np.exp(1j*phi)]])
def wp_r(phi,th):
    return np.dot(rot(-th),np.dot(wp(phi),rot(th)))
def pol():
    return np.array([[1,0],[0,0]])
def pol_r(th):
    return np.dot(rot(-th),np.dot(pol(),rot(th)))

j0 = np.array([1,0])
retardance = 2*pi/4.0

pol_angle  = np.linspace(0,2*pi,100)

cnt = 0

phi = np.linspace(0,pi,500)
vis = phi*0

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
    
plt.plot(phi/pi,vis)

plt.xlim([0,1])
plt.grid(True)
plt.ylim([0,1.1])

plt.show()