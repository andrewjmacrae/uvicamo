import numpy as np
from matplotlib import pyplot as plt


f_path = 'data/Nov14_2019_data/'
t = np.loadtxt(f_path+'time_dat.csv')*1000
idx = np.loadtxt(f_path+'index_pulse.csv')
dat = np.loadtxt(f_path+'pol_data.csv')

fig,(ax1,ax2) = plt.subplots(1,2,figsize=[12,4])
ax1.plot(t,idx,label = 'Trigger Channel')
ax1.plot(t,dat,label = 'PD data')

def extract_triggers(trig_dat,thrsh=1):
    trigz = np.array([])
    for d in range(len(trig_dat)-1):
        if trig_dat[d+1]-trig_dat[d] > thrsh:
            trigz = np.append(trigz,int(d))
    return trigz.astype(int)

trigz = extract_triggers(idx)
ax1.plot(t[trigz+1],idx[trigz+1],'o',label='Detected Triggers')
dt = t[1]-t[0]

for k in range(len(trigz)-1):
	chunk = dat[trigz[k]:trigz[k+1]]
	L = len(chunk)
	tt = np.linspace(0,dt*L,L)
	print(str(k)+'th chunck had '+str(L)+' points')
	ax2.plot(tt,chunk)

ax1.legend(loc='best')
ax1.set_xlabel('Time [ms]')
ax2.set_xlabel('Time [ms]')
ax1.set_ylabel('Pulse Detection')
ax2.set_ylabel('PD signal')
ax1.set_title('Trigger Events')
ax2.set_title('Extracted Chunks')
plt.show()