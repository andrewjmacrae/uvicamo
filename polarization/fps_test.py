from matplotlib import pyplot as plt, animation
import numpy as np
import time

# globals
cnt = 0
Npts = int(1e2)
t = np.linspace(0,1,Npts)

# Set up plots
fig,ax  = plt.subplots(1,1,figsize = [5,5])
ln, = ax.plot(np.random.random(Npts))
then = time.time()
fps = 30


print('Requesting '+str(fps)+' fps ...')

def init():
	ax.set_xlim([min(t),max(t)])
	return ln,

# Dummy function to return dummy data
def get_dat_data(enn):
	return np.random.random(enn)

# Runs at a rate of fps frames per second.
# Not sure what happens if get_dat_data takes longer than a frame ...
def animate(idx):
	global cnt,then,fps
	cnt+=1
	ln.set_data(t,get_dat_data(Npts))
	if cnt%fps == 0: # Update count 
		now = time.time()
		dt = now-then
		then = now
		print(str(round(fps/dt))+' FPS')
	return ln,

ani = animation.FuncAnimation(fig, animate,init_func = init, interval=1000/fps)
plt.show()