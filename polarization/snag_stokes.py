import sys
import csv
import swptools as swp
import os.path
import json
import numpy as np

#f_base = '/'
f_base = '/home/pi/Documents/data/Dec12/'

csv_name = 'Fresnel_stokes.csv'
f_name = f_base+csv_name

if len(sys.argv) == 2:
    deg = int(sys.argv[1])
    print('Writing to file: '+f_name)
else:
    print('Wrong number of args, b')
    sys.exit()
    
# This code snags a trace

from daqhats import mcc118, OptionFlags, HatIDs, HatError
from daqhats_utils import select_hat_device, enum_mask_to_string, chan_list_to_mask

sim_settings_file = 'settings/simsettings.json'
daq_settings_file = 'settings/daqsettings.json'
swp_settings_file = 'settings/swpsettings.json'

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

hat.a_in_scan_start(channel_mask, samples_per_channel, scan_rate, options)
read_result = hat.a_in_scan_read(samples_per_channel, timeout)

y1 = np.array(read_result.data[::2]) - bg_level
y2 = read_result.data[1::2]  

hat.a_in_scan_stop()
hat.a_in_scan_cleanup()

# Next need to convert the CHONK to a set o' Stokes veX

trigz = swp.extract_triggers(y2)
Nchunks = len(trigz)-1

S = np.zeros(4)

for k in range(Nchunks):
    chunk = np.array(y1[trigz[k]:trigz[k+1]])
    S += swp.get_stokes_from_chunk(chunk,wp_ret = wp_phi,phs_ofst = trigger_phase,verbose = False)
S/=Nchunks
S /= S[0]
DOP = np.sqrt(S[1]**2 + S[2]**2 + S[3]**2)

print(f'DOP = {DOP}')

with open(f_name, 'a', newline = '\n') as cfile:
    wrtr = csv.writer(cfile,delimiter = ',')
    wrtr.writerow([deg,S[0],S[1],S[2],S[3]])

