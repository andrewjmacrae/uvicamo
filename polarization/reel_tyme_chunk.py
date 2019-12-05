# First task is to get the reading and plotting multiple channels -DONE (indexing is weird)
# Then determine time for reading (without plotting) to get bound on fps - DONE (less than 5 ms lag)
# Then get continuous scan DONE (had to use stop() and clear() or it barfed)
# Then Apply chunking alogorithm (Done)
# Then extract Fourier coefficients via integration (Done)
# Then calibrate phase !! (Done)
# Then extract Stokes Params and plot Pol. ellipse (Done)
# Then plot bar graph (Done) and poincare (Not so Done)

# minimum number of triggers
# minimum amount of intensity

from daqhats import mcc118, OptionFlags, HatIDs, HatError
from daqhats_utils import select_hat_device, enum_mask_to_string, chan_list_to_mask
import numpy as np
from matplotlib import pyplot as plt, animation

phs = .05

timeout=.5
channels = [0, 2]
channel_mask = chan_list_to_mask(channels)
num_channels = len(channels)

address = select_hat_device(HatIDs.MCC_118)
hat = mcc118(address)

samples_per_channel = 1000
scan_rate = 5000.0
# options = OptionFlags.DEFAULT
options = OptionFlags.CONTINUOUS

address = select_hat_device(HatIDs.MCC_118)
hat = mcc118(address)

Nsmp = samples_per_channel
Ts = 1/scan_rate
tTot = Ts*Nsmp
t = np.linspace(0,tTot-Ts,Nsmp)

#create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))

ln1, = ax1.plot([], [], lw=2)

#initialize bar graph as a variable outside animation function so
#  we don't have to clear out data each frame
bar = plt.bar([0, 1, 2, 3], [0, 0, 0, 0], align='center')

def polarization_ellipse(S0,S1,S2):
    '''
    given a stokes vector, this function plots the corresponding
    polarization ellipse

    returns a list of x and y values on the Ex-Ey plane
    
    NOTE: S0 must be equal to 1 and S1,S2 should be normalized by S0
    '''

    #solve for psi, the angle from the x axis of the ellipse
    
    '''
    old method before arctan2 was implemented
    if S1 == 0:
        psi = np.pi/4
    else:
        psi = 0.5*np.arctan2(S2,S1)
    '''
    
    psi = 0.5*np.arctan2(S2,S1)

    #define ellipse parameters from stokes vectors
    a = np.sqrt(0.5*(S0+np.sqrt(S1**2+S2**2)))
    b = np.sqrt(0.5*(S0-np.sqrt(S1**2+S2**2)))
    rot = np.matrix([[np.cos(psi), -1*np.sin(psi)],
                     [np.sin(psi), np.cos(psi)]])
    ba = b/a
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
    ax1.set_xlabel('Ex')
    ax1.set_ylabel('Ey')
    ax1.set_title('Polarization Ellipse')
    ax1.grid(True)
    ax1.set_aspect('equal')

    #graph parameters for Bar graph
    ax2.set_xlim(-0.6, 3.6)
    ax2.set_ylim(-1.05, 1.05)
    #ax2.set_ylim(-10, 10)
    plt.xticks([0, 1, 2, 3], ['S0', 'S1', 'S2', 'S3'])
    ax2.set_xlabel('Stokes Parameter')
    ax2.set_ylabel('Value of Stokes Parameter')
    ax2.set_title('Stokes Parameters')
    
    print('Phase: '+str(round(phs*180/2*np.pi,1))+' deg')
    return ln1,
    
def extract_triggers(trig_dat,thrsh=1):
    trigz = np.array([])
    for d in range(len(trig_dat)-1):
        if trig_dat[d+1] - trig_dat[d] > thrsh:
            trigz = np.append(trigz,int(d))
    return trigz.astype(int)

def animate_fun(idx):
    global phs
    
    hat.a_in_scan_start(channel_mask, samples_per_channel, scan_rate, options)
    read_result = hat.a_in_scan_read(samples_per_channel, timeout)
    y1 = read_result.data[::2]
    y2 = read_result.data[1::2]
    
    trigz = extract_triggers(y2)
    
    print('found '+str(len(trigz))+' chunks')
    
    hat.a_in_scan_stop()
    hat.a_in_scan_cleanup()
    
    C2w,S2w,C4w,S4w,s0,C = 0,0,0,0,0,0
    Nchunks = len(trigz)
    
    print('y1:',y1)
    print('mean(y1)',np.mean(y1))
    
    for k in range(Nchunks-1):
        y = y1[trigz[k]:trigz[k+1]]        
        wt = np.linspace(0,2*np.pi,len(y))
        C2w += np.trapz(y*np.cos(2*wt + phs))
        S2w += np.trapz(y*np.sin(2*wt + phs))
        C4w += np.trapz(y*np.cos(4*wt + phs))
        S4w += np.trapz(y*np.sin(4*wt + phs))
        s0 += np.mean(y)
        C+=np.trapz(y)
    
    #computing Stokes Parameters from Fourier Shenanigans
    #S0 = C
    #S0 = s0
    S1 = 2*C4w
    S2 = 2*S4w
    S3 = S2w
    S0 = np.trapz(y1)*2 - S1/2
    
    if S0**2<S1**2+S2**2+S3**2:
        print('DO NOT TRUST THE FOLLOWING RESULTS')
        
    if np.mean(y1) < 0.01:
        #just by eye for now
        #this value is not independent from gain of the detector
        #how about S/N? This would require a noise reading before the experiment
            #or would it... ?
        print('NOT ENOUGH LIGHT FOR ACCURATE RESULTS')
        
    if Nchunks < 3:
        #Accurate results can actually be acquired with as few as 2 chunks
        #at 2 chunk limit small rpm fluctuations can lead to nan's
        #accuracy increases as nchunks increases, there may be a better way
            #ie %diff from convergent value? 
        print('WAVEPLATE IS NOT ROTATING FAST ENOUGH FOR ACCURATE DATA')
        
    print('nchunks', Nchunks)
    
    #should normalize via S0 but giving funky results...
    nrm = np.sqrt(S1**2 + S2**2 + S3**2)
    #nrm  = S0
    S0/=nrm
    S1/=nrm
    S2/=nrm
    S3/=nrm
  
    print('S = ['+str(round(S0,3))+','+str(round(S1,3))+','+str(round(S2,3))+','+str(round(S3,3))+']')
    
    print('nrm',nrm)
    
    # degree of polarization

    print("DOP:", np.sqrt(abs(S1**2+S2**2+S3**3)))
    
    #S = [1,S1,S2,S3]
    S = [S0,S1,S2,S3]
    x,y = polarization_ellipse(1,S1,S2)
    
    
    ln1.set_data(x,y)
    
    for i in range(len(S)):
        bar[i].set_height(S[i])
    
    return ln1, bar,
    
annie = animation.FuncAnimation(fig,animate_fun,init_func = init_animation, interval = 250)
plt.show()