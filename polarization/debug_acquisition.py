# This tests the acquisition, processing, and stokes extraction algorithms
# We ain't got no time for no fancy animation here. 


import numpy as np
from matplotlib import pyplot as plt, animation


we_live_in_a_simulation = True
sim = we_live_in_a_simulation
pi = np.pi

if not sim:
    from daqhats import mcc118, OptionFlags, HatIDs, HatError
    from daqhats_utils import select_hat_device, enum_mask_to_string, chan_list_to_mask

phs = .05

samples_per_channel = 1000
scan_rate = 5000.0
    
if not sim:
    timeout=.5
    channels = [0, 2]
    channel_mask = chan_list_to_mask(channels)
    num_channels = len(channels)
    address = select_hat_device(HatIDs.MCC_118)
    hat = mcc118(address)

    # options = OptionFlags.DEFAULT
    options = OptionFlags.CONTINUOUS

    address = select_hat_device(HatIDs.MCC_118)
    hat = mcc118(address)

Nsmp = samples_per_channel
Ts = 1/scan_rate
tTot = Ts*Nsmp
t = np.linspace(0,tTot-Ts,Nsmp)

# In this program, the left plot shows the entire plot plus triggers, the right shows the chunked data
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

ln_raw, = ax1.plot([], [], lw=2)
ln_trg, = ax1.plot([], [], lw=2)
ln_cnk, = ax2.plot([], [], lw=2)
ln_wtf, = ax2.plot([], [], '--',lw=2)


def init_animation():    
    
    #graph parameters for Ellipse
    ax1.set_xlim(0, tTot*1000)
    ax1.set_ylim(-0.1, 5.1)
    ax1.set_xlabel('Time [ms]')
    ax1.set_ylabel('Signal Level [V]')
    ax1.set_title('Raw Data')
    ax1.grid(True)

    #graph parameters for Bar graph
    ax2.set_xlim(0, 1)
    ax2.set_ylim(-.1, 2.1)
    
    ax2.set_xlabel('Theta $\\theta/2\\pi$')
    ax2.set_ylabel('Signal Level')
    ax2.set_title('ChunX')
    
    return ln_raw,ln_trg,ln_cnk,ln_wtf,

def sim_pol_data(S0,w0,t0,sig_level=1,ns_level = 0):
    global phs
    Npts = len(t0)
    a = (2*S0[0]+S0[1])/4
    b = S0[3]/2
    c = S0[1]/4
    d = S0[2]/4
    ns = ns_level*np.random.randn(Npts)
    return a + b*np.sin(2*w0*t0 + phs) + c*np.cos(4*w0*t0 + phs) + d*np.sin(4*w0*t0 + phs) + ns


def extract_triggers(trig_dat,thrsh=1):
    trigz = np.array([])
    for d in range(len(trig_dat)-1):
        if trig_dat[d+1] - trig_dat[d] > thrsh:
            trigz = np.append(trigz,int(d)+1)
    return trigz.astype(int)

def animate_fun(idx):
    global phs,t
    if not sim:
        hat.a_in_scan_start(channel_mask, samples_per_channel, scan_rate, options)
        read_result = hat.a_in_scan_read(samples_per_channel, timeout)
        y1 = read_result.data[::2]
        y2 = read_result.data[1::2]
    else:
        DP = 1
        w = 2*np.pi*3500/60
        S = 2*np.array([1,0,0,1])
        y1 = sim_pol_data(S,w,t,ns_level=.01)
        y2 = 5*(np.mod(w*t,2*np.pi) < np.pi/12)
    
    trigz = extract_triggers(y2)
    
    
    if not sim:
        hat.a_in_scan_stop()
        hat.a_in_scan_cleanup()

    Nchunks = len(trigz)

    # C0w, C2w, S2w, C4w, S4w, Mnw = 0,0,0,0,0,0
    # print(f'Working with {Nchunks} chunks')
    # # print('y1:',y1)
    # # print('mean(y1)',np.mean(y1))
    
    # for k in range(Nchunks-1):
    #     y = y1[trigz[k]:trigz[k+1]]
    #     wt = np.linspace(0,2*np.pi,len(y))
    #     C2w += np.trapz(y*np.cos(2*wt + phs, wt))
    #     S2w += np.trapz(y*np.sin(2*wt + phs, wt))
    #     C4w += np.trapz(y*np.cos(4*wt + phs, wt))
    #     S4w += np.trapz(y*np.sin(4*wt + phs, wt))
    #     Mnw += np.mean(y)
    #     C0w += np.trapz(y,wt)

    #     print(f'Yo! C2 = {np.trapz(y*np.cos(2*wt + phs, wt))} and it should be zero!')

    k0 = int(np.random.rand()*(Nchunks-1))
    y0 = y1[trigz[k0]:trigz[k0+1]]
    t0 = np.linspace(0,2*pi,len(y0))

    ln_raw.set_data(t*1e3,y1)
    ln_trg.set_data(t*1e3,y2)
    ln_cnk.set_data(t0/(2*pi),y0)

    wtf = np.linspace(0,2*np.pi,len(y0))
    ln_wtf.set_data(wtf/(2*pi),(1+np.sin(2*wtf + phs)))
    
    print(f'Chunk {k0}, cosine: {round(np.trapz(np.cos(2*wtf + phs)*y0,wtf)/pi,2)}, sine: {round(np.trapz(np.sin(2*wtf + phs)*y0,wtf)/pi,2)}')

    return ln_raw,ln_trg,ln_cnk,ln_wtf,
    
annie = animation.FuncAnimation(fig,animate_fun,init_func = init_animation, interval = 250)
plt.show()
