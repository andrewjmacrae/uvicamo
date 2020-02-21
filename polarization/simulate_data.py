import numpy as np
from matplotlib import animation, pyplot as plt

pi = np.pi

fig,ax = plt.subplots()

DP = .8

def sim_pol_data(S0,w0,t0,sig_level=1,ns_level = 0):
    Npts = len(t)
    a = (2*S0[0]+S0[1])/4
    b = S0[3]/2
    c = S0[1]/4
    d = S0[2]/4
    ns = ns_level*np.random.randn(Npts)
    return a + b*np.sin(2*w*t) + c*np.cos(4*w*t) + d*np.sin(4*w*t) + ns

w = 2*pi*88.8
t = np.linspace(0,2*pi/w*8.2,5000)


S = 2*np.array([1,.706*DP,0*DP,-.706*DP])
y = sim_pol_data(S,w,t)
triggered = 5*(np.mod(w*t,2*pi) < pi/12)
ln1, = ax.plot(t,y)
ln2, = ax.plot(t,triggered)
ax.set_xlim(min(t),max(t))

def init_animation():
	return ln1,ln2,

def animate(frm):
	DP = .8
	S = 2*np.array([1,DP*np.cos(frm/18.)/np.sqrt(2),DP*np.sin(frm/18.)/np.sqrt(2),-.706*DP])

	ln1.set_data(t,sim_pol_data(S,w,t,ns_level=.05))
	ln2.set_data(t,5*(np.mod(w*t,2*pi) < pi/12))
	return ln1,ln2,
ani = animation.FuncAnimation(fig,animate,init_func = init_animation, interval = 125)

plt.show()