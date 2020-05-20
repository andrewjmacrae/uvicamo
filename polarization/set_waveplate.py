import json
import numpy as np

# Enter max, min, and bg voltage here
v_max = 2.84
v_min = 1.74
v_bg = 0.04

eta = (v_min-v_bg)/(v_max-v_bg)

with open('settings/swpsettings.json','r') as f:
	params = json.load(f)
	print(f'Waveplate phase set from '+str(round(params['wp_phi'],3)) + ' to ' + str(round(eta,3)))
	params['wp_phi'] = eta

with open('settings/swpsettings.json','w') as f:
	json.dump(params, f)