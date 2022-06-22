import numpy as np

def get_stokes_from_chunk(cnk,wp_ret = np.pi/2,phs_ofst = 0,verbose = False):
    """
    Gets stokes vectors from a given chunk of data
    
    Parameters:
    	cnk: numpy array
    		chunk data value, represents full rotation of waveplate
    		
    	wp_ret: float (optional)
    		retardance of waveplate. Default np.pi/2
    		
    	phs_ofset: float (optional)
    		offset angle of the waveplate WRT the horizontal. Default 0
    		
    	verbose: bool (optional)
    		enables error messaging. Default False
    		
    Returns:
    	numpy array, stokes vector
    """
    
    a0,b0,c0,d0,n0 = 0,0,0,0,0

    wt = np.linspace(0,2*np.pi,len(cnk))
    a0 = np.trapz(cnk,wt)/(2*np.pi)
    n0 = np.trapz(cnk*np.cos(2*(wt-phs_ofst)),wt)/np.pi
    b0 = np.trapz(cnk*np.sin(2*(wt-phs_ofst)),wt)/np.pi
    c0 = np.trapz(cnk*np.cos(4*(wt-phs_ofst)),wt)/np.pi
    d0 = np.trapz(cnk*np.sin(4*(wt-phs_ofst)),wt)/np.pi

    cd = np.cos(wp_ret)
    sd = np.sin(wp_ret)

    S0 = 2*(a0 - c0*(1+cd)/(1-cd))
    S1 = 4*c0/(1-cd)
    S2 = 4*d0/(1-cd)
    S3 = -2*b0/sd
    
    nrm = S0
    
    if nrm == 0:
        if verbose:
            print('Error! S0 = 0. Something went terribly wrong')
        nrm=1 #This is probably sketchy
    
    if n0 > np.sqrt(S1**2 + S2**2 + S3**2)*1e-3 and verbose:
        print(f'Warning, large sin(2w) conponent detected ({n0}). Check alignment!')
    
    return np.array([S0,S1,S2,S3])/nrm


def extract_triggers(trig_dat,thrsh=1,schmidt = 10):
    """
    Gives an array of timestamps of each rotation trigger
    
    Parameters:
    	trig_dat: numpy array of integers
    		raw data from the trigger channel
    	
    	thrsh: float (optional)
    		voltage threshold a signal must pass to be considered a good trigger. Default 1
    	
    	schmidt: int (optional)
    		uses schmidt triggers to avoid multiple triggerings caused by noise. Default 10
    		
    Returns:
    	array of integers, timestamps of rotation triggers
    """
    
    trigz = np.array([])
    deadzone = 0
    for d in range(len(trig_dat)-1):
        deadzone = max(0,deadzone-1)
        if trig_dat[d+1] - trig_dat[d] > thrsh and deadzone is 0:
            trigz = np.append(trigz,int(d))
            deadzone = schmidt
    return trigz.astype(int)

def get_polarization_ellipse(S,n_points = 200,scale_by_dop = True, verbose = False):
    """
    From a stokes vector, finds points to plot the resulting polarization ellipse
    
    Parameters:
    	S: numpy array
    		stokes vector
    	
    	n_points: int (optional)
    		number of points in x,y arrays, equal to number of points in ellipse. Default 200
    	
    	scale_by_dop: bool (optional)
    		scales down semimajor axis of elipse in accordance with degree of polarization. Default True
    		# make that one easier to understand, or just shorten, if possible, (Scale down ellipse for partial polarization)?
    	
    	verbose: bool (optional)
    		turns on error messages. Default False
    		
    Returns:
    	x,y numpy arrays of coordiantes for graphing the ellipse
    """
    
    S/=S[0]
    
    DPol = np.sqrt(sum(S[1:]**2))
    for k in range(len(S)):
        if abs(S[k]) > 1:
            if verbose:
                print(f'WARNING: S[{k}] is greater than 100% Clipping to 1.')
            S[k] = S[k]/abs(S[k])

    if scale_by_dop:
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
    x = np.linspace(-a, a, n_points)

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

def sim_pol_data(S0,w0,t0,sig_level=1,ns_level = 0,digitize_mV = 0,v_bias = 0,dphi = np.pi/2,ofst = 0):
    """
    Simulates what would have been recorded by the polarimeter, based on a mock data json.
    Allows testing the polarimeter without inputing real data
    
    Parameters:
    	S0: numpy array
    		initial stokes vector to base data on [rephrase]
    	
    	w0: float
    		frequency of the waveplate
    	
    	t0: float
    		current time
    	
    	sig_level: float (optional)
    		max amplitude of the signal. Default 1
    	
    	ns_level: float (optional)
    		amount of noise simulated. Default 0
    	
    	digitize_mV: float (optional)
    		digitization level. Default 0
    	
    	v_bias: float (optional)
    		DC voltage offset. Default 0
    	
    	dphi: float (optional)
    		retardance of the waveplate. Default np.pi/2
    	
    	ofst: float (optional)
    		phase offset. Default 0
    		
    Returns:
    	trc: 
    """
    
    Npts = len(t0)
    a = S0[0]/2 + (1+np.cos(dphi))*S0[1]/4
    b = -S0[3]*np.sin(dphi)/2
    c = (1-np.cos(dphi))*S0[1]/4
    d = (1-np.cos(dphi))*S0[2]/4
    ns = ns_level*np.random.randn(Npts)
    trc =  (a + b*np.sin(2*w0*t0 - 4*np.pi*ofst) + c*np.cos(4*w0*t0 - 8*np.pi*ofst) + d*np.sin(4*w0*t0 - 8*np.pi*ofst))*sig_level + ns + v_bias
    if digitize_mV > 0:
        trc = np.around(trc*1000/digitize_mV)*digitize_mV/1000
    return trc
 
