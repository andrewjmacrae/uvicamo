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

def extract_triggers(trig_dat,thrsh=1):
    trigz = np.array([])
    for d in range(len(trig_dat)-1):
        if trig_dat[d+1] - trig_dat[d] > thrsh:
            trigz = np.append(trigz,int(d))
    return trigz.astype(int)

w = 2*pi*88.8
t = np.linspace(0,2*pi/w*8.2,5000)


S = 2*np.array([1,.706*DP,0*DP,-.706*DP])
y = sim_pol_data(S,w,t)
triggered = 5*(np.mod(w*t,2*pi) < pi/12)

# for k in range(len(trigz)-1):
#     chonk = y[trigz[k]:trigz[k+1]]
#     print(len(chonk))
#     tt = np.linspace(0,2*pi,len(chonk))
#     print(f'(C,S) = ({round(np.trapz(chonk*np.cos(2*tt),tt),3)}, {round(np.trapz(chonk*np.sin(2*tt),tt),3)})')
ln1, = ax.plot(t,y)
ln2, = ax.plot(t,triggered)
ax.set_xlim(min(t),max(t))

def init_animation():
	return ln1,ln2,

def animate(frm):
    DP = .8
    S = 2*np.array([1,DP*np.cos(frm/18.)/np.sqrt(2),DP*np.sin(frm/18.)/np.sqrt(2),-.706*DP])
    y1 = sim_pol_data(S,w,t,ns_level=.05)
    y2 = 5*(np.mod(w*t,2*pi) < pi/12)
    trigz = extract_triggers(y2)
    for k in range(len(trigz)-1):
        chonk = y[trigz[k]:trigz[k+1]]
        print(len(chonk))
        tt = np.linspace(0,2*pi,len(chonk))
        print(f'(C,S) = ({round(np.trapz(chonk*np.cos(2*tt),tt),3)}, {round(np.trapz(chonk*np.sin(2*tt),tt),3)})')
    trigz = extract_triggers(y2)
    ln1.set_data(t,y1)
    ln2.set_data(t,y2)
    return ln1,ln2,

ani = animation.FuncAnimation(fig,animate,init_func = init_animation, interval = 125)

plt.show()