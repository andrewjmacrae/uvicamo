# First task is to get the reading and plotting multiple channels -DONE (indexing is weird)
# Then determine time for reading (without plotting) to get bound on fps - DONE (less than 5 ms lag)
# Then get continuous scan DONE (had to use stop() and clear() or it barfed)
# Then Apply chunking alogorithm (Done)
# Then extract Fourier coefficients via integration (Done)
# Then calibrate phase !! (Done)
# Then extract Stokes Params and plot Pol. ellipse (Done)
# Then plot bar graph (Done) and poincare (Not so Done)
# add raw values to plots (done)
# scale ellipse by degree of polarization (done)

#clean up
#DOP in a function
#scale S123 by S0
#make sure raw values refresh on plot (done)


# minimum number of triggers
# minimum amount of intensity

# Add simulation Mode, debug algorithm for (simulated) known system (done)

we_live_in_a_simulation = True
sim = we_live_in_a_simulation

poincare = True
trace = poincare


if not sim:
    from daqhats import mcc118, OptionFlags, HatIDs, HatError
    from daqhats_utils import select_hat_device, enum_mask_to_string, chan_list_to_mask
import numpy as np
from matplotlib import pyplot as plt, animation
if poincare:
    from mpl_toolkits.mplot3d import Axes3D

phs = .05*(1-sim)

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

#create figure
if trace:
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    ax3.axis('off')
    ax3 = fig.add_subplot(1,3,3, projection = '3d')
    ln3, = ax3.plot([],[], lw=2, label = 'trace')
else:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5),\
                                   sharey = True,\
                                   gridspec_kw={'wspace': 0.1})

ln1, = ax1.plot([], [], lw=2)

#text variables to update
txt1 = ax1.text(-.95,0.9,'',fontsize = 12)
# txt2 = ax1.text(-1, 0.85, '', fontsize = 12)

#initialize bar graph as a variable outside animation function so
#  we don't have to clear out data each frame
bar = ax2.bar([0, 1, 2, 3], [0, 0, 0, 0], align='center')

def polarization_ellipse(S0,S1,S2,DPol):
    '''
    given a stokes vector, this function plots the corresponding
    polarization ellipse

    returns a list of x and y values on the Ex-Ey plane
    
    NOTE: S0 must be equal to 1 and S1,S2 should be normalized by S0
    '''

    #solve for psi, the angle from the x axis of the ellipse
    
    psi = 0.5*np.arctan2(S2,S1)
    DPol = 1
    #define ellipse parameters from stokes vectors
    a = np.sqrt(0.5*(S0+np.sqrt(S1**2+S2**2)))*DPol
    b = np.sqrt(0.5*(S0-np.sqrt(S1**2+S2**2)))*DPol
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
    x = x1+x2[::-1]
    y = y1+y2[::-1]

    return x, y

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
        
        # define points for poincare sphere wire frame
        u = np.linspace(0.0,2*np.pi,90)
        v = np.linspace(0.0,np.pi,180)
        lu = np.size(u)
        lv = np.size(v)
        x = np.zeros((lu,lv))
        y = np.zeros((lu,lv))
        z = np.zeros((lu,lv))
        for uu in range(0,lu):
            for vv in range(0,lv):
                x[uu,vv]= np.cos(u[uu])*np.sin(v[vv])
                y[uu,vv]= np.sin(u[uu])*np.sin(v[vv])
                z[uu,vv]= np.cos(v[vv])
    
        persist=15 #no of fading out points to include
        trans=np.arange(1.0,0.0,-1.0/(persist+1))

        ax3.plot_wireframe(x, y, z, rstride=9, cstride=15, linewidth=3, color='gray', alpha=0.25)
        
        ax3.axis('off')
        ax3.set_xlim(-0.6,0.6)
        ax3.set_ylim(-0.6,0.6)
        ax3.set_zlim(-0.6,0.6)
        ax3.set_title('Poincaré')
        ax3.set_aspect('equal')
        
        return ln1,bar,txt1,ln3,
    
    
    print('Phase: '+str(round(phs*180/2*np.pi,1))+' deg')
    return ln1,bar,txt1,

def sim_pol_data(S0,w0,t0,sig_level=1,ns_level = 0):
    Npts = len(t0)
    a = (2*S0[0]+S0[1])/4
    b = S0[3]/2
    c = S0[1]/4
    d = S0[2]/4
    ns = ns_level*np.random.randn(Npts)
    return a + b*np.sin(2*w0*t0) + c*np.cos(4*w0*t0) + d*np.sin(4*w0*t0) + ns

def extract_triggers(trig_dat,thrsh=1):
    trigz = np.array([])
    for d in range(len(trig_dat)-1):
        if trig_dat[d+1] - trig_dat[d] > thrsh:
            trigz = np.append(trigz,int(d))
    return trigz.astype(int)

def animate_fun(idx):
    global phs,t
    if not sim:
        hat.a_in_scan_start(channel_mask, samples_per_channel, scan_rate, options)
        read_result = hat.a_in_scan_read(samples_per_channel, timeout)
        y1 = read_result.data[::2]
        y2 = read_result.data[1::2]
    else:
        DP = .65
        Phi = float(idx/18.)
        w = 2*np.pi*1400/60
        S = 3*np.array([1,DP*np.cos(Phi)/np.sqrt(2),DP*np.sin(Phi)/np.sqrt(2),DP/np.sqrt(2)])
        y1 = sim_pol_data(S,w,t,ns_level=.00)
        y2 = 5*(np.mod(w*t,2*np.pi) < np.pi/12)
    
    trigz = extract_triggers(y2)
        
    if not sim:
        hat.a_in_scan_stop()
        hat.a_in_scan_cleanup()

# I'm rewriting this part of the algorithm, if only because I'm too stoopid to git it. - AM
    a0,n0,b0,c0,d0 = 0,0,0,0,0
    Nchunks = len(trigz)-1

    Nroll = 0 # Holds rolling average of pts per chunk (PPC)
    
    for k in range(Nchunks):
        chunk = y1[trigz[k]:trigz[k+1]]
        Nroll = (Nroll*k + len(chunk))/(k+1) # Update (PPC)
        wt = np.linspace(0,2*np.pi,len(chunk))
        a0 += np.trapz(chunk,wt)/(2*np.pi*Nchunks)
        n0 += np.trapz(chunk*np.cos(2*wt),wt)/(np.pi*Nchunks)
        b0 += np.trapz(chunk*np.sin(2*wt),wt)/(np.pi*Nchunks)
        c0 += np.trapz(chunk*np.cos(4*wt),wt)/(np.pi*Nchunks)
        d0 += np.trapz(chunk*np.sin(4*wt),wt)/(np.pi*Nchunks)

    #computing Stokes Parameters from Fourier Shenanigans
    S0 = 2*(a0-c0)
    S1 = 4*c0
    S2 = 4*d0
    S3 = 2*b0

    nrm = S0
    S = np.array([S0,S1,S2,S3])/nrm
    
    DOP = np.sqrt(S[1]**2 + S[2]**2 + S[3]**2)

# All kinds of warnings!!!
    if Nroll < 180:
        print('Warning: insufficient points per revolution for accurate data - slow\'er down!')
    if DOP-1 > .03:
        print(f'Warning: Possible unphysical DOP = {round(DOP,2)} measured...')
    if n0/nrm > 1e-2:
        print(f'Warning: Possible alignment error: large cos(2wt) component detected (S,C) = ({b0/nrm},{n0/nrm})')
        
    if np.mean(y1) < 0.08:
        #just by eye for now
        #this value is not independent from gain of the detector
        #how about S/N? This would require a noise reading before the experiment
            #or would it... ?
        print('Warning: Low light level detected ...')
        
    if Nchunks < 3:
        #Accurate results can actually be acquired with as few as 2 chunks
        #at 2 chunk limit small rpm fluctuations can lead to nan's
        #accuracy increases as nchunks increases, there may be a better way
            #ie %diff from convergent value? 
        print('Warning: Insufficient periods. Is waveplate spinning?')
      
    # print(f'S = {np.around(S,3)}')
        
    x,y = polarization_ellipse(S[0],S[1],S[2],DOP)
    
    txt1.set_text(f'DOP: {round(DOP,3)}')
    # txt2.set_text(f'Mean Signal: {np.mean(y1)}')
        
    ln1.set_data(x,y)
    
    for i in range(len(S)):
        bar[i].set_height(S[i])

    if trace:
        #ln3.set_data(range(len(y1)),y1)
        ax3.scatter3D(S[1],S[2],S[3], alpha = 0.4, color = 'b')
        
        return ln1, bar, txt1, ln3,
    
    return ln1, bar, txt1,
    
annie = animation.FuncAnimation(fig,animate_fun,init_func = init_animation, interval = 250)
plt.show()
