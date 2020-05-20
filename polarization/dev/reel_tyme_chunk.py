# reel_tyme_chunk.py
# A. MacRae, A. McKay, S. Wilkenson

we_live_in_a_simulation = True
sim = we_live_in_a_simulation

if not sim:
    from daqhats import mcc118, OptionFlags, HatIDs, HatError
    from daqhats_utils import select_hat_device, enum_mask_to_string, chan_list_to_mask
import numpy as np
from matplotlib import pyplot as plt, animation

# ------- Simulation data ----
sim_digitize = 1000*20/(2**12) # MCC118 is 12bit and +/- 10 V
sim_siglevel = 1
sim_ns_level = 0.03
sim_DOP = .7
sim_vbias = 0.00
wp_phase = 1.982 # percentage deviation from perfect QWP
sim_phase_offset = 0.404
S_sim = np.array([1,0,0,1])
#----

trace_debug_mode = True
trace = trace_debug_mode

phs = .404
wp_phi = np.arccos(-.4)

samples_per_channel = 1000
scan_rate = 20000.0
auto_scale_y_trace = False
    
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

#create figure
if trace:
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    ln3, = ax3.plot([],[], lw=2, label = 'trace')
else:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5),\
                                   sharey = True,\
                                   gridspec_kw={'wspace': 0.1})

ln1, = ax1.plot([], [], lw=2)

#text variables to update
txt1 = ax1.text(-.95,0.9,'',fontsize = 12)
txt2 = ax1.text(-.95, 0.8, '', fontsize = 12, color = 'blue')
txt_err = ax1.text(-1.25,-1.35,'', fontsize = 10, color = 'red')
#initialize bar graph as a variable outside animation function so
#  we don't have to clear out data each frame
bar = ax2.bar([0, 1, 2, 3], [0, 0, 0, 0], align='center')
bar2 = ax2.bar([0, 1, 2, 3], [0, 0, 0, 0], align='edge',alpha = .3)

def polarization_ellipse(S):
    '''
    given a stokes vector, this function plots the corresponding
    polarization ellipse

    returns an array of x and y values on the Ex-Ey plane
    '''    
    
    S/=S[0]
    
    DPol = np.sqrt(sum(S[1:]**2))
    for k in range(len(S)):
        if abs(S[k]) > 1:
            S[k] = S[k]/abs(S[k])

    S1 = S[1]/DPol
    S2 = S[2]/DPol
    S3 = S[3]/DPol

    psi = 0.5*np.arctan2(S2,S1)
    chi = 0.5*np.arcsin(S3)
    
    a = 1
    b = np.tan(chi)
    
    ba = b/a
    rot = np.matrix([[np.cos(psi), -1*np.sin(psi)],
                     [np.sin(psi), np.cos(psi)]])
    x1, x2, y1, y2 = [], [], [], []
    #create an x array for plotting ellipse y values
    x = np.linspace(-a, a, 200)

    for x in x:
        #cartesian equation of an ellipse
        Y1 = ba*np.sqrt(a**2-x**2)
        #Y1 reflection about the x-axis
        #rotate the ellipse by psi
        XY1 = np.matrix([[x],
                         [Y1]])
        XY2 = np.matrix([[x],
                         [-Y1]])
        y1.append(float((rot*XY1)[1]))
        x1.append(float((rot*XY1)[0]))
        y2.append(float((rot*XY2)[1]))
        x2.append(float((rot*XY2)[0]))

    #x2,y2 reversed in order so that there is continuity in the ellipse (no line through the middle)
    x = (x1+x2[::-1])
    y = (y1+y2[::-1])

    return np.array(x)*DPol, np.array(y)*DPol

def init_animation():    
    
    #graph parameters for Ellipse
    ax1.set_xlim(-1.0, 1.0)
    ax1.set_ylim(-1.0, 1.0)
    ax1.set_xlabel('$E_x$')
    ax1.set_ylabel('$E_y$')
    ax1.set_title('Polarization Ellipse')
    ax1.grid()
    ax1.set_aspect('equal')
    
    #add unit circle around ellipse
    t_circ = np.linspace(0,2*np.pi,1000)
    x = [np.cos(T) for T in t_circ]
    y = [np.sin(T) for T in t_circ]
    ax1.plot(x,y, color = 'gray', linestyle = '--', alpha = 0.5)

    #graph parameters for Bar graph
    ax2.set_xlim(-0.6, 3.6)
    ax2.set_ylim(-1.05, 1.05)
    ax2.set_xticks([0, 1, 2, 3])
    ax2.set_xticklabels(['S0', 'S1', 'S2', 'S3'])
    ax2.set_xlabel('Stokes Parameter')
    #ax2.set_ylabel('Value of Stokes Parameter')
    ax2.set_title('Stokes Parameters')
    
    if trace:
        ax3.set_xlim(-1,1000)
        ax3.set_ylim(0,3)
        ax3.set_title('Trace')
        ax3.grid()
        
        return ln1,bar,txt1,ln3,txt_err
    
    print('Phase: '+str(round(phs*180/2*np.pi,1))+' deg')
    return ln1,bar,txt1,

def sim_pol_data(S0,w0,t0,sig_level=1,ns_level = 0,digitize_mV = 0,v_bias = 0,dphi = np.pi/2,ofst = 0):
    Npts = len(t0)
    a = S0[0]/2 + (1+np.cos(dphi))*S0[1]/4
    b = -S0[3]*np.sin(dphi)/2
    c = (1-np.cos(dphi))*S0[1]/4
    d = (1-np.cos(dphi))*S0[2]/4
    ns = ns_level*np.random.randn(Npts)
    trc =  (a + b*np.sin(2*w0*t0 + ofst) + c*np.cos(4*w0*t0 + ofst) + d*np.sin(4*w0*t0 + ofst))*sig_level + ns + v_bias
    if digitize_mV > 0:
        trc = np.around(trc*1000/digitize_mV)*digitize_mV/1000
    return trc

def extract_triggers(trig_dat,thrsh=1):
    trigz = np.array([])
    for d in range(len(trig_dat)-1):
        if trig_dat[d+1] - trig_dat[d] > thrsh:
            trigz = np.append(trigz,int(d))
    return trigz.astype(int)

def animate_fun(idx):
    global phs,t, S_sim
    estr = 'warnings: '
    if sim:
        DP = sim_DOP
        Phi = float(idx/18.)
        w = 2*np.pi*5100/60        
        if np.mod(int(idx/50),3) == 0:
            S_sim = 3*np.array([1,DP*np.cos(Phi)/np.sqrt(2),DP*np.sin(Phi)/np.sqrt(2),DP/np.sqrt(2)])
            estr = 'ellip-pol. '+estr
        elif np.mod(int(idx/50),3) == 1:            
            S_sim = 3*np.array([1,DP*np.cos(Phi),DP*np.sin(Phi),0])
            estr = 'lin-pol. '+estr
        else:
            S_sim = 3*np.array([1,0,0,DP])
            estr = 'circ-pol. '+estr
        y1 = sim_pol_data(S_sim,w,t,ns_level=sim_ns_level,sig_level = sim_siglevel,digitize_mV=sim_digitize, v_bias = sim_vbias,dphi=wp_phase,ofst = sim_phase_offset)
        y2 = 5*(np.mod(w*t,2*np.pi) < np.pi/12)
    else:
        hat.a_in_scan_start(channel_mask, samples_per_channel, scan_rate, options)
        read_result = hat.a_in_scan_read(samples_per_channel, timeout)
        y1 = read_result.data[::2]
        y2 = read_result.data[1::2]  

    trigz = extract_triggers(y2)
    
    if not sim:
        hat.a_in_scan_stop()
        hat.a_in_scan_cleanup()

# I'm rewriting this part of the algorithm, if only because I'm too stoopid to git it. - AM
    a0,n0,b0,c0,d0 = 0,0,0,0,0
    Nchunks = len(trigz)-1

    Nroll = 0 # Holds rolling average of pts per chunk (PPC). Used for accuracy warning.
    
    # For each chunk we calculate the 0,2w, and 4w comonents, averaged over all chunks
    for k in range(Nchunks):
        chunk = y1[trigz[k]:trigz[k+1]]
        Nroll = (Nroll*k + len(chunk))/(k+1) # Update (PPC)
        wt = np.linspace(0,2*np.pi,len(chunk))
        a0 += np.trapz(chunk,wt)/2
        n0 += np.trapz(chunk*np.cos(2*wt+phs),wt)
        b0 += np.trapz(chunk*np.sin(2*wt+phs),wt)
        c0 += np.trapz(chunk*np.cos(4*wt+phs),wt)
        d0 += np.trapz(chunk*np.sin(4*wt+phs),wt)

    #computing Stokes Parameters from Fourier Shenanigans (see analysis document)
    cd = np.cos(wp_phi)
    sd = np.sin(wp_phi)
    prf = 1./(Nchunks*np.pi)

    S0 = 2*prf*(a0 - c0*(1+cd)/(1-cd))
    S1 = 4*prf*c0/(1-cd)
    S2 = 4*prf*d0/(1-cd)
    S3 = -2*prf*b0/sd

    nrm = S0

    n0*=prf/nrm
    S = np.array([S0,S1,S2,S3])/nrm
    
    DOP = np.sqrt(S[1]**2 + S[2]**2 + S[3]**2)

    # estr = 'warnings: '
# All kinds of warnings!!!
    if Nroll < 180:
        # print('Warning: insufficient points per revolution for accurate data - slow\'er down!')
        estr+=f'PPC too low ({int(Nroll)})    '
    if DOP-1 > .03:
        # print(f'Warning: Possible unphysical DOP = {round(DOP,2)} measured...')
        estr += f'Unphysical DOP ({round(DOP,3)})    '
    if n0 > 4e-2:
        # print(f'Warning: Possible alignment error: large cos(2wt) component detected (S,C) = ({b0/nrm},{n0/nrm})')
        estr += f'non-zero cos(2w): {round(n0,2)} cf {round(b0*prf/nrm,2)}    '
        
    if np.mean(y1) < 0.08:
        #just by eye for now
        #this value is not independent from gain of the detector
        #how about S/N? This would require a noise reading before the experiment
            #or would it... ?
        # print('Warning: Low light level detected ...')
        estr += f'Light level too low    '
        
    if Nchunks < 3:
        #Accurate results can actually be acquired with as few as 2 chunks
        #at 2 chunk limit small rpm fluctuations can lead to nan's
        #accuracy increases as nchunks increases, there may be a better way
            #ie %diff from convergent value? 
        # print('Warning: Insufficient periods. Is waveplate spinning?')
        estr += f'Insufficient chunks    '
      
    # print(f'S = {np.around(S,3)}')
    
    ###################DEBUG ZONE#####################
    
    # print(f'Avg PPC: {Nroll}')
    
    ###################################################
    txt_err.set_text(estr)
    x,y = polarization_ellipse(S)
    
    txt1.set_text(f'DOP: {round(DOP,3)}')
    if sim:
        urr = 100*abs(DOP - sim_DOP)/sim_DOP
        txt2.set_text(f'Error: {round(urr,1)}%')
    else:
        txt2.set_text(f'Mean Signal: {np.mean(y1)}')
        
    ln1.set_data(x,y)
    
    for i in range(len(S)):
        bar[i].set_height(S[i])
        if sim:
            bar2[i].set_height(S_sim[i]/S_sim[0])

    if trace:
        ln3.set_data(range(len(y1)),y1)
        ln3.set_data(range(len(chunk)),chunk)
        if auto_scale_y_trace:
            ax3.set_ylim(min(chunk) + 0.001, max(chunk) +0.001)
        ax3.set_xlim(0, len(chunk))
        
        return ln1, bar, txt1, ln3,
    
    return ln1, bar, txt1,
    
annie = animation.FuncAnimation(fig,animate_fun,init_func = init_animation, interval = 250)
plt.show()
