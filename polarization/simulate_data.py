import numpy as np
from matplotlib import animation, pyplot as plt

fig,ax, = plt.subplots()

DP = .8

def sim_pol_data(S0,w0,t0,sig_level=1,ns_level = 0):
    Npts = len(t)
    a = (2*S0[0]+S0[1])/4
    b = S0[3]/2
    c = S0[1]/4
    d = S0[2]/4
    ns = 0.01*np.random.randn(Npts)
    return a + b*np.sin(2*w*t) + c*np.cos(4*w*t) + d*np.sin(4*w*t) + ns

w = 2*np.pi*88.8
t = np.linspace(0,2*np.pi/w*8.2,5000)


S = 2*np.array([1,.706*DP,0*DP,-.706*DP])
y = sim_pol_data(S,w,t)

ax.plot(y)

plt.show()