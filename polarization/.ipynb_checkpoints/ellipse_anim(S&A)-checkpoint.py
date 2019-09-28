### https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import cmath
import time
# From stokes parameters, looking to create the ellipse and then animate it on a raspPi
#to then be used for the GUI once we start getting real time # data.
# start = time.time()
def polarization_ellipse(S):
    #set individual stokes parameters from stokes vector
    S0,S1,S2,S3 = S[0],S[1],S[2],S[3]

    #solve for psi, the angle from the x axis of the ellipse
    if S1 == 0:
        psi = np.pi/4
    else:
        psi = 0.5*np.arctan(S2/S1)

    #define ellipse parameters from stokes vectors
    a=np.sqrt(0.5*(S0+np.sqrt(S1**2+S2**2)))
    b=np.sqrt(0.5*(S0-np.sqrt(S1**2+S2**2)))
    ba = b/a
    rot = np.matrix([[np.cos(psi),-1*np.sin(psi)],
                     [np.sin(psi),np.cos(psi)]])
    x1,x2,y1,y2 = [],[],[],[]
    #create an x array for plotting ellipse y values
    x = np.linspace(-a, a, 200)
#     print(a)
 
    for x in x: 
        Y1 = ba*np.sqrt(a**2-x**2)
        #relection about the x-axis
        Y2 = -Y1
#         Y2 = -ba*np.sqrt(a**2-x**2)
        #rotate the ellipse by psi
        XY1 = np.matrix([[x],
                        [Y1]])
        XY2 = np.matrix([[x],
                        [Y2]])
        y1.append(float((rot*XY1)[1]))
        x1.append(float((rot*XY1)[0]))
        y2.append(float((rot*XY2)[1]))
        x2.append(float((rot*XY2)[0]))
    
    #x2,y2 reversed in order so that there is continuity in the ellipse (no line through the middle)
    x=x1+x2[::-1]
    y=y1+y2[::-1]

    return x,y


fig = plt.figure(figsize = (10,10))
ax = plt.axes(xlim=(-1.0, 1.0), ylim=(-1.0, 1.0))
ax.set_xlabel('Ex')
ax.set_ylabel('Ey')
ax.set_title('Polarization Ellipse')
plt.grid()
line, = ax.plot([], [], lw=2)

def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
#later, iterate through S instead of phi
def animate(phi): 
    xanim = []
    yanim = []
    S = [1,0,np.cos(phi),np.sin(phi)]
    xanim,yanim = polarization_ellipse(S)
    line.set_data(xanim,yanim)
    
    return line,


# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=np.linspace(1, 2*np.pi, 800), interval=1, blit=True)
    
plt.show()
