from daqhats import mcc118, OptionFlags, HatIDs, HatError
from daqhats_utils import select_hat_device, enum_mask_to_string, chan_list_to_mask
import numpy as np
from matplotlib import pyplot as plt

#address = select_hat_device(HatIDs.MCC_118)
address = 0
hat = mcc118(address)

read_request_size = 500
timeout = 5.0



channels = [0]
channel_mask = chan_list_to_mask(channels)
num_channels = len(channels)

samples_per_channel = 10000
scan_rate = 100000.0 # Breaks above 100,000
options = OptionFlags.DEFAULT
 
 
actual_scan_rate = hat.a_in_scan_actual_rate(num_channels, scan_rate)
hat.a_in_scan_start(channel_mask, samples_per_channel, scan_rate,options)

read_result = hat.a_in_scan_read(read_request_size, timeout)

print(len(read_result.data))

y = read_result.data

N = len(y)

Ts = 1/scan_rate

t = np.linspace(0,(N-1)*Ts,N)
f = np.linspace(0,scan_rate,N)

Y = abs(np.fft.fft(y-np.mean(y))*2/N)
idx = np.argmax(Y[0:int(N/2)])



fig,(ax1,ax2) = plt.subplots(1,2,figsize=[12,5])
ax1.plot(1e6*t,y,'o-')
ax1.set_xlabel('Time [$\mu$s]')
ax1.set_ylabel('Voltage [V]')
ax1.set_xlim([0,1000])
ax1.set_title('Input Sine Wave: A = 1.7V, f = 15.0 kHz')
ax2.plot(f/1000,Y,label = 'FFT')
ax2.plot(f[idx]/1000,Y[idx],'o',label = 'max at '+str(round(f[idx]/1000,3))+' kHz')
ax2.grid(True)
ax2.legend()
ax2.set_title('Frequency Spectrum (Fs = 100 kHz)')
ax2.set_xlim([0,scan_rate/2000])
ax2.set_xlabel('Frequency [kHz]')
ax2.set_ylabel('Spectral Power')
plt.show()