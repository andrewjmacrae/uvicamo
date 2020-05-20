#gen_default_json.py
import os.path
import json

sim_settings_file = 'settings/simsettings.json'
daq_settings_file = 'settings/daqsettings.json'
swp_settings_file = 'settings/swpsettings.json'


# ------- Simulation data ----
if os.path.isfile(sim_settings_file):
	print(f'Found simulation data file: {sim_settings_file}... Skipping')
else:
	sim_dict = {	'sim_digitize': 1000*20/(2**12),
					'sim_siglevel': 1,
					'sim_ns_level': 0.03,
					'sim_DOP': 0.7,
					'sim_vbias': 0.003,
					'wp_phase': 1.982,
					'sim_phase_offset': 0.404,
					'sim_poltype': 'right'
				}
	print(f'No simulation file found. Creating file: {sim_settings_file}')
	with open(sim_settings_file,'w') as f:
		json.dump(sim_dict, f)

# ------- SWP Settings Data ----
if os.path.isfile(swp_settings_file):
	print(f'Found SWP data file: {swp_settings_file}... Skipping')
else:
	swp_dict = {	'trigger_phase': 0.404,
					'wp_phi': 1.982,
					'bg_level': 0.003,
					'auto_scale_y_trace': False,
					'log_data_file': '',
					'protect_overwrite': True
				}
	print(f'No SWP file found. Creating file: {swp_settings_file}')
	with open(swp_settings_file,'w') as f:
		json.dump(swp_dict, f)

# ------- DAQ Settings Data ----
if os.path.isfile(daq_settings_file):
	print(f'Found DAQ data file: {daq_settings_file}... Skipping')
else:
	daq_dict = {	'samples_per_channel': 1000,
					'scan_rate': 20000,
					'channels': [0,2],
                                        'timeout': 5
				}
	print(f'No DAQ json file found. Creating file: {daq_settings_file}')
	with open(daq_settings_file,'w') as f:
		json.dump(daq_dict, f)

print('Finished')
