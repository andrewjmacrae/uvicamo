# ps_polvis.py

run_offline = False

import swptools as swp
import numpy as np
from matplotlib import pyplot as plt, animation
from mpl_toolkits.mplot3d import axes3d
import json
import os.path

if not run_offline:
    from daqhats import mcc118, OptionFlags, HatIDs, HatError
    from daqhats_utils import select_hat_device, enum_mask_to_string, chan_list_to_mask	

sim_settings_file = 'settings/simsettings.json'
daq_settings_file = 'settings/daqsettings.json'
swp_settings_file = 'settings/swpsettings.json'

do_save = False
# Load data from json file
if run_offline:
	if not os.path.isfile(sim_settings_file):
		print(f'Error: simulation file {sim_settings_file} not found.')
		print('Run \'gen_default_json.py\' to generate default file first.')
		exit()
	else:
		with open(sim_settings_file,'r') as f:
			simpams = json.load(f)
			sim_digitize = simpams['sim_digitize']
			sim_siglevel = simpams['sim_siglevel']
			sim_ns_level = simpams['sim_ns_level']
			sim_DOP = simpams['sim_DOP']
			sim_bg_level = simpams['sim_vbias']
			sim_wp_phi = simpams['wp_phase']
			sim_trigger_phase = simpams['sim_phase_offset']
			sim_poltype = simpams['sim_poltype']
			if sim_poltype == 'right':
				S_sim = np.array([1,0,0,1])
			elif sim_poltype == 'lin':
				S_sim = np.array([1,1,0,0])
			else:
				S_sim = np.array([1,np.sqrt(.3),np.sqrt(.3),np.sqrt(.4)])

if not os.path.isfile(daq_settings_file):
	print(f'Error: simulation file {sim_settings_file} not found.')
	print('Run \'gen_default_json.py\' to generate default file first.')
	exit()
else:
	with open(daq_settings_file,'r') as f:
		daqpams = json.load(f)
		samples_per_channel = daqpams['samples_per_channel']
		scan_rate = daqpams['scan_rate']
		channels = daqpams['channels']
		timeout = daqpams['timeout']
		Nsmp = samples_per_channel
		Ts = 1/scan_rate
		tTot = Ts*Nsmp
		t = np.linspace(0,tTot-Ts,Nsmp)
		if not run_offline:
			channel_mask = chan_list_to_mask(channels)
			num_channels = len(channels)
			address = select_hat_device(HatIDs.MCC_118)
			hat = mcc118(address)
			options = OptionFlags.CONTINUOUS
            
if not os.path.isfile(swp_settings_file):
	print(f'Error: settings file {swp_settings_file} not found.')
	print('Run \'gen_default_json.py\' to generate default file first.')
	exit()
else:
	with open(swp_settings_file,'r') as f:
		swppams = json.load(f)
		trigger_phase = swppams['trigger_phase']
		wp_phi = swppams['wp_phi']
		auto_scale_y_trace = swppams['auto_scale_y_trace']
		bg_level = swppams['bg_level']
		data_log_file = swppams['log_data_file']
		if data_log_file != '':
			do_save = True
			with open(data_log_file,'w') as f:
				f.write('')

# Set up plot canvas
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
# Plots
ln1, = ax1.plot([], [], lw=2)
bar = ax2.bar([0, 1, 2, 3], [0, 0, 0, 0], align='center')
bar2 = ax2.bar([0, 1, 2, 3], [0, 0, 0, 0], align='edge',alpha = .3)
#`ln3, = ax3.plot([],[], lw=2, label = 'trace')
ax3 = fig.add_subplot(1,3,3, projection = '3d')
#text variables to update
txt1 = ax1.text(-.95,0.9,'',fontsize = 12)
txt2 = ax1.text(-.95, 0.8, '', fontsize = 12, color = 'blue')
txt_err = ax1.text(-1.25,-1.35,'', fontsize = 10, color = 'red')

pt = ax3.scatter([0.5],[0.5],[0.707],facecolor='tab:blue',s=100)
ln3, = ax3.plot(np.array([0,0.5]),np.array([0,0.0]),np.array([0,0.0]),color='tab:red',lw=3)


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
    # Now we begin the munumental task of creating the 3D plot:
    # Begin by defining a set of points for the sphere.
    u = np.linspace(0.0,2*np.pi,30)
    v = np.linspace(0.0,np.pi,60)
    lu = np.size(u)
    lv = np.size(v)
    X = np.zeros((lu,lv))
    Y = np.zeros((lu,lv))
    Z = np.zeros((lu,lv))

    for uu in range(0,lu):
        for vv in range(0,lv):
            X[uu,vv]= np.cos(u[uu])*np.sin(v[vv])
            Y[uu,vv]= np.sin(u[uu])*np.sin(v[vv])
            Z[uu,vv]= np.cos(v[vv])
    srf = ax3.plot_surface(X, Y, Z,alpha=.2,color = 'gray')
    # Now add the axes and spherical gridlines:
    tht = np.linspace(0,2*np.pi,360)
    ax3.plot([-1,1],[0,0],[0,0],lw=2,color = 'black')
    ax3.plot([0,0],[-1,1],[0,0],lw=2,color = 'black')
    ax3.plot([0,0],[0,0],[-1,1],lw=2,color = 'black')
    ax3.plot(np.cos(tht),np.sin(tht),0*tht,color = 'gray',alpha = 0.6)
    
    ax3.plot(np.cos(tht),0*tht,np.sin(tht),color = 'black',alpha = 0.5)
    ax3.plot(0*tht,np.sin(tht),np.cos(tht),color = 'black',alpha = 0.5)
    ax3.plot(np.sin(tht),np.cos(tht),0*tht,color = 'black',alpha = 0.5)

    # Add the labels
    fsz = 16
    ax3.text(1.3,0,0,'H',color='tab:blue',fontsize=fsz)
    ax3.text(-1.3,0,0,'V',color='tab:blue',fontsize=fsz)
    ax3.text(0,-1.3,0,'-45',color='tab:red',fontsize=fsz)
    ax3.text(0,1.1,0,'+45',color='tab:red',fontsize=fsz)
    ax3.text(0,0,1.2,'R',color='tab:green',fontsize=fsz)
    ax3.text(0,0,-1.3,'L',color='tab:green',fontsize=fsz)

    ax3.view_init(elev=30, azim=30)

    for phi0 in np.linspace(-np.pi/2,np.pi/2,12):
        ax3.plot(np.cos(tht)*np.cos(phi0),np.sin(tht)*np.cos(phi0),np.sin(phi0),color = 'gray',alpha = 0.4) 
    return ln1,bar,txt1,txt_err
    
def animate_fun(idx):
    global phs,t, S_sim
    estr = 'warnings: '
    if run_offline:
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
        y1 = swp.sim_pol_data(S_sim,w,t,ns_level=sim_ns_level,sig_level = sim_siglevel,digitize_mV=sim_digitize, v_bias = sim_bg_level,dphi=sim_wp_phi,ofst = sim_trigger_phase)
        y2 = 5*(np.mod(w*t,2*np.pi) < np.pi/12)
    else:
        hat.a_in_scan_start(channel_mask, samples_per_channel, scan_rate, options)
        read_result = hat.a_in_scan_read(samples_per_channel, timeout)
        y1 = np.array(read_result.data[::2]) - bg_level
        y2 = read_result.data[1::2]
        hat.a_in_scan_stop()
        hat.a_in_scan_cleanup()

    trigz = swp.extract_triggers(y2)

    Nchunks = len(trigz)-1

    Nroll = 0 # Holds rolling average of pts per chunk (PPC). Used for accuracy warning.
    # For each chunk we calculate the 0,2w, and 4w comonents, averaged over all chunks
    # get_stokes_from_chunk(cnk,wp_ret = np.pi/2,phs_ofst = 0,verbose = False):

    S = np.zeros(4)

    for k in range(Nchunks):
        chunk = np.array(y1[trigz[k]:trigz[k+1]])
        #get_stokes_from_chunk(cnk,wp_ret = np.pi/2,phs_ofst = 0,verbose = False):
        S += swp.get_stokes_from_chunk(chunk,wp_ret = wp_phi,phs_ofst = trigger_phase,verbose = False)
        Nroll = (Nroll*k + len(chunk))/(k+1) # Update (PPC)
    #computing Stokes Parameters from Fourier Shenanigans (see analysis document)
    S/=Nchunks
    S /= S[0]
    
    DOP = np.sqrt(S[1]**2 + S[2]**2 + S[3]**2)
    if do_save:
	    with open(data_log_file,'a') as f:
	    	f.write(f'{S[1]},{S[2]},{S[3]},{DOP}\n')

    # estr = 'warnings: '
# All kinds of warnings!!!
    if Nroll < 180:
        # print('Warning: insufficient points per revolution for accurate data - slow\'er down!')
        estr+=f'PPC too low ({int(Nroll)})    '
    if DOP-1 > .03:
        # print(f'Warning: Possible unphysical DOP = {round(DOP,2)} measured...')
        estr += f'Unphysical DOP ({round(DOP,3)})    '
    # if n0 > 4e-2:
    #     # print(f'Warning: Possible alignment error: large cos(2wt) component detected (S,C) = ({b0/nrm},{n0/nrm})')
    #     estr += f'non-zero cos(2w): {round(n0,2)} cf {round(b0*prf/nrm,2)}    '
        
    if np.mean(y1) < 0.08:
        estr += f'Light level too low    '
        
    if Nchunks < 3:
        estr += f'Insufficient chunks    '
      
    txt_err.set_text(f'({round(np.mean(y1),2)},{round(S[1],2)},{round(S[2],2)},{round(S[3],2)})\n'+estr)
    x,y = swp.get_polarization_ellipse(S)
    
    txt1.set_text(f'DOP: {round(DOP,3)}')
    if run_offline:
        urr = 100*abs(DOP - sim_DOP)/sim_DOP
        txt2.set_text(f'Error: {round(urr,1)}%')
    else:
        txt2.set_text(f'Mean Signal: {np.mean(y1)}')
        
    ln1.set_data(x,y)
    
    for i in range(len(S)):
        bar[i].set_height(S[i])
        if run_offline:
            bar2[i].set_height(S_sim[i]/S_sim[0])

    D = np.sqrt(S[1]**2 + S[2]**2)
    if D < 1e-3:
        D = 1e-7
    else:
        D = np.sqrt(1-S[3]**2)/D


    x0 = S[1]*D
    y0 = S[2]*D
    z0 = S[3]
    
    ln3.set_data(np.array([0,x0]),np.array([0,y0]))
    ln3.set_3d_properties(np.array([0,z0]))
    pt._offsets3d = ([x0],[y0],[z0])
    return ln1, bar, txt1, 
    
annie = animation.FuncAnimation(fig,animate_fun,init_func = init_animation, interval = 250)
plt.show()
