# Visualization of polarization rotation through sugar water FTW!!!
import numpy as np
from matplotlib import animation, pyplot as plt

pi = np.pi


fig, (ax, ax2) = plt.subplots(1,2,figsize=(12, 4))
ax.set_aspect('equal')

frq = 4
n_frames = 200
print('Using '+str(n_frames)+' frames.')


ax.set_xlim(-1.2,1.2)
ax.set_ylim(-1.2,1.2)
ax.set_label('Phase = 0')
ax.grid(True)


lX = 0.5
rX = 0.75
x = np.linspace(lX,rX,10)
ax2.fill_between(x,-.5,.5,alpha=.5,color = '#ff7f0e')

ax2.set_xlim(0,1)
ax2.set_ylim(-.5,.5)


qax = ax.quiver([0,0,0],[0,0,0],[1,1,1],[1,1,1],angles='xy', scale_units='xy', scale=1.)
# ray = ax2.quiver([],)
line, = ax2.plot([0,1],[0,0],lw=2)

# lines = [line]

def init():
#     lines[0].set_data(0,0)
#     qax = ax.quiver(0,0,1,1)
    return qax,ax,ax2,

def animate(i):
    idx = 0
    phi0 = 0
    if i > n_frames/2:
        idx = i-n_frames/2
    if i > 3*n_frames/4:
        idx = n_frames/4 # 3N/4 - N/4
    tht = frq*i*2*pi/n_frames
    phs = phi0+.75*idx*2*pi/n_frames
    ExR = np.cos(tht)/2
    ExL = np.cos(-tht-phs)/2
    ExTot = ExR+ExL
    EyR = np.sin(tht)/2
    EyL = np.sin(-tht-phs)/2
    EyTot = EyR+EyL
    
    qax.set_UVC([ExR,ExL,ExTot],[EyR,EyL,EyTot],[100,100,0])
    line.set_data([0,i*1.0/n_frames],[0,0]) 

    ax.set_title('Polarization angle:'+str(round(phs*90/pi))+' deg')

    return qax,ax,ax2,

anim = animation.FuncAnimation(fig,animate,init_func=init,frames=n_frames,interval=30,blit=False)

plt.show()