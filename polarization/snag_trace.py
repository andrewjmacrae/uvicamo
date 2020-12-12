# run_cal_phase.py

# Grab a trace with several periods and determine phase offset by minimizing cos(2wt) term

import numpy as np
from matplotlib import pyplot as plt
import csv

f_base = '/home/pi/Documents/data/Dec12/'
csv_name = 'fresnel_trc.csv'

from daqhats import mcc118, OptionFlags, HatIDs, HatError
from daqhats_utils import select_hat_device, enum_mask_to_string, chan_list_to_mask

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

hat.a_in_scan_stop()
hat.a_in_scan_cleanup()


with open(f_base+csv_name, 'w', newline = '\n') as cfile:
    wrtr = csv.writer(cfile,delimiter = ',')
    for k in range(Nsmp):
        wrtr.writerow([t[k],y1[k],y2[k]])
print('wrote file: '+f_base+csv_name)

fig, ax = plt.subplots(1,1,figsize = [7,4])
ax.plot(t,y1,label= 'Trace')
ax.plot(t,y2,label='Triggers')
ax.legend()
ax.grid(True)
plt.show()
