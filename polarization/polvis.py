# ps_polvis.py

run_offline = True # testing (True) or measuring (False)
do_save = False

import swptools as swp
import numpy as np
from matplotlib import pyplot as plt, animation
from mpl_toolkits.mplot3d import axes3d
import json
import os.path
if not run_offline:
    from daqhats import mcc118, OptionFlags, HatIDs, HatError
    from daqhats_utils import select_hat_device, enum_mask_to_string, chan_list_to_mask	

# applying variable names to the paths to the json files
sim_settings_file = 'settings/simsettings.json' # settings for running offline (simulation)
daq_settings_file = 'settings/daqsettings.json' # settings for daqhats_utils
swp_settings_file = 'settings/swpsettings.json' # settings for swptools


# load sim settings json
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


# load daq settings json
if not os.path.isfile(daq_settings_file):
	print(f'Error: simulation file {daq_settings_file} not found.')
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


# load swp settings json
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
# Define plots
ln1, = ax1.plot([], [], lw=2)
bar = ax2.bar([0, 1, 2, 3], [0, 0, 0, 0], align='center')
bar2 = ax2.bar([0, 1, 2, 3], [0, 0, 0, 0], align='edge',alpha = .3)
ln3, = ax3.plot([],[], lw=2, label = 'trace')
# Define text variables to update
txt1 = ax1.text(-.95,0.9,'',fontsize = 12)
txt2 = ax1.text(-.95, 0.8, '', fontsize = 12, color = 'blue')
txt_err = ax1.text(-1.25,-1.35,'', fontsize = 10, color = 'red')

def init_animation():
	"""
	Sets the parameters for the 3 graphs, used as the initial function in FuncAnimation.
	Allows setting the graph parameters without cluttering up a single function

	Parameters:
		None

	Returns:
		ln1, bar, ln3: paramaterized functions
		txt1, txt_err: accompanying text
	"""
	
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
	ax2.set_title('Stokes Parameters')
	
	# graph parameters for trace
	ax3.set_xlim(-1,1000)
	ax3.set_ylim(0,3)
	ax3.set_title('Trace')
	ax3.grid()
	
	return ln1,bar,txt1,ln3,txt_err

def animate_fun(idx):
	"""
	Defines the function FuncAnimation will animate.
	
	Parameters:
		idx: int, current frame in the animation
		# note: this parameter is not to be called independantly, its inclusion is
			needed because FuncAnimation requires such a parameter to function

	Returns:
		ln1, bar, ln3: 
		graphs to be animated

		txt1: 
	
	"""
	
	global phs, t, S_sim
	
	# setup error message text
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
			
		y1 = swp.sim_pol_data(S_sim, w, t,
				      ns_level = sim_ns_level,
				      sig_level = sim_siglevel,
				      digitize_mV = sim_digitize, 
				      v_bias = sim_bg_level,
				      dphi = sim_wp_phi,
				      ofst = sim_trigger_phase)
		
		y2 = 5*(np.mod(w*t,2*np.pi) < np.pi/12)
		
	else:
		hat.a_in_scan_start(channel_mask, samples_per_channel, scan_rate, options)
		read_result = hat.a_in_scan_read(samples_per_channel, timeout)
		y1 = np.array(read_result.data[::2]) - bg_level
		y2 = read_result.data[1::2]
		hat.a_in_scan_stop()
		hat.a_in_scan_cleanup()

	trigz = swp.extract_triggers(y2)

	# I'm rewriting this part of the algorithm, if only because I'm too stoopid to git it. - AM
	Nchunks = len(trigz)-1

	Nroll = 0 # Holds rolling average of pts per chunk (PPC). Used for accuracy warning.
	# For each chunk we calculate the 0,2w, and 4w comonents, averaged over all chunks

	S = np.zeros(4)

	for k in range(Nchunks):
		chunk = np.array(y1[trigz[k]:trigz[k+1]])
		S += swp.get_stokes_from_chunk(chunk,wp_ret = wp_phi,phs_ofst = trigger_phase,verbose = False)
		Nroll = (Nroll*k + len(chunk))/(k+1) # Update (PPC)
		
	#computing Stokes Parameters from Fourier Shenanigans (see analysis document)
	S/=Nchunks
	S /= S[0]

	DOP = np.sqrt(S[1]**2 + S[2]**2 + S[3]**2)
	if do_save:
		with open(data_log_file,'a') as f:
			f.write(f'{S[1]},{S[2]},{S[3]},{DOP}\n')
	
	
	# checking for errors, and appending error messages to estr as needed
	if Nroll < 180:
		estr+=f'PPC too low ({int(Nroll)})    '
	if DOP-1 > .03:
		estr += f'Unphysical DOP ({round(DOP,3)})    '

	if np.mean(y1) < 0.08:
		estr += f'Light level too low    '

	if Nchunks < 3:
		estr += f'Insufficient chunks    '
	
	
	# setting up text
	txt_err.set_text(f'({round(np.mean(y1),2)},{round(S[1],2)},{round(S[2],2)},{round(S[3],2)})\n'+estr)
	txt1.set_text(f'DOP: {round(DOP,3)}')
	
	if run_offline:
		urr = 100*abs(DOP - sim_DOP)/sim_DOP
		txt2.set_text(f'Error: {round(urr,1)}%')
	else:
		txt2.set_text(f'Mean Signal: {np.mean(y1)}')
	
	
	# application of data to graphs
	x,y = swp.get_polarization_ellipse(S)
	
	ln1.set_data(x,y)

	for i in range(len(S)):
		bar[i].set_height(S[i])
		if run_offline:
			bar2[i].set_height(S_sim[i]/S_sim[0])

	ln3.set_data(range(len(y1)),y1)
	ax3.set_xlim(0, len(y1))
	
	
	return ln1, bar, txt1, ln3,
    
annie = animation.FuncAnimation(fig,animate_fun,init_func = init_animation, interval = 250)
plt.show()
