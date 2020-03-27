run_cal_phase.py

# Grab a trace with several periods and determine phase offset by minimizing cos(2wt) term

import numpy as np
from matplotlib import pyplot as plt
from daqhats import mcc118, OptionFlags, HatIDs, HatError
from daqhats_utils import select_hat_device, enum_mask_to_string, chan_list_to_mask

def extract_triggers(trig_dat,thrsh=1):
    trigz = np.array([])
    for d in range(len(trig_dat)-1):
        if trig_dat[d+1] - trig_dat[d] > thrsh:
            trigz = np.append(trigz,int(d))
    return trigz.astype(int)

samples_per_channel = 1000
scan_rate = 20000.0

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

hat.a_in_scan_start(channel_mask, samples_per_channel, scan_rate, options)
read_result = hat.a_in_scan_read(samples_per_channel, timeout)

y1 = read_result.data[::2]
y2 = read_result.data[1::2]  

trigz = extract_triggers(y2)

hat.a_in_scan_stop()
hat.a_in_scan_cleanup()

enns = array([])

phases = linspace(0,2*pi,1000)

for phs in phases:
	for k in range(Nchunks):
        chunk = y1[trigz[k]:trigz[k+1]]
        wt = np.linspace(0,2*np.pi,len(chunk))
        n0 += np.trapz(chunk*np.cos(2*wt+phs),wt)
    enns = append(enns,n0/Nchunks)

fig, (ax1,ax2) = subplots(1,2,figsize = [13,4])

ax1.plot(phases,enns)
ax2.plot(t,y1,lw=2)
ax2.plot(t,y2,'--',lw=1)

plt.show()