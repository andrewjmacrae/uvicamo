### https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import cmath

#TO DO:
#
#grab stokes parameters from output text file in real time

# From stokes parameters, looking to create the ellipse and then animate it on a raspPi
#to then be used for the GUI once we start getting real time # data.


def polarization_ellipse(S):
    #set individual stokes parameters from stokes vector
    S0, S1, S2, S3 = S[0], S[1], S[2], S[3]

    #solve for psi, the angle from the x axis of the ellipse
    if S1 == 0:
        psi = np.pi/4
    else:
        psi = 0.5*np.arctan(S2/S1)

    #define ellipse parameters from stokes vectors
    a = np.sqrt(0.5*(S0+np.sqrt(S1**2+S2**2)))
    b = np.sqrt(0.5*(S0-np.sqrt(S1**2+S2**2)))

    x1, x2, y1, y2 = [], [], [], []
    #create an x array for plotting ellipse y values
    x = np.linspace(-5, 5, 1500)

    for x in x:
        if a**2-x**2 >= 0:  # to avoid negative values in the square root
            #cartesian equation of an ellipse
            Y1 = (b/a)*np.sqrt(a**2-x**2)
            #relection about the x-axis
            Y2 = -1*(b/a)*np.sqrt(a**2-x**2)
            #rotate the ellipse by psi
            XY1 = np.matrix([[x],
                             [Y1]])
            XY2 = np.matrix([[x],
                             [Y2]])
            rot = np.matrix([[np.cos(psi), -1*np.sin(psi)],
                             [np.sin(psi), np.cos(psi)]])

            y1.append(float((rot*XY1)[1]))
            x1.append(float((rot*XY1)[0]))
            y2.append(float((rot*XY2)[1]))
            x2.append(float((rot*XY2)[0]))

    #x2,y2 reversed in order so that there is continuity in the ellipse (no line through the middle)
    x = x1+x2[::-1]
    y = y1+y2[::-1]

    return x, y


#create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))

#set graph parameters for Ellipse
ax1.set_xlim(-1.0, 1.0)
ax1.set_ylim(-1.0, 1.0)
ax1.set_xlabel('Ex')
ax1.set_ylabel('Ey')
ax1.set_title('Polarization Ellipse')
ax1.grid(True)
ax1.set_aspect('equal')

#set graph parameters for Stokes bar graph
ax2.set_xlim(-0.6, 3.6)
ax2.set_ylim(-1.05, 1.05)
plt.xticks([0, 1, 2, 3], ['S0', 'S1', 'S2', 'S3'])
ax2.set_xlabel('Stokes Parameter')
ax2.set_ylabel('Value of Stokes Parameter')
ax2.set_title('Stokes Parameters')

line1, = ax1.plot([], [], lw=2)


def init():
    line1.set_data([], [])

    return line1,


#initialize bar graph as a global variable, outside animation function so we don't have to clear out data each frame
bar = plt.bar([0, 1, 2, 3], [0, 0, 0, 0], align='center')

# animation function.  This is called sequentially
#later, iterate through S instead of phi


def animate(phi):
    x1, y1, x2, y2 = [], [], [], []
    S = [1, 0, np.cos(phi), np.sin(phi)]

    #replace above line with grabbing from text file

    #Method 1: read in all and grab last value
    stokes = np.genfromtxt(stokes_output.txt, dtype=float)
    S = S[-1]

    #Method 2: Not sure what function can do this but read in last line
    #  such that it is removed from text file

    x1, y1 = polarization_ellipse(S)
    x2 = [0, 1, 2, 3]
    y2 = S

    line1.set_data(x1, y1)
    #bar = plt.bar(x2,y2, align = 'center')
    for i in range(len(S)):
        bar[i].set_height(S[i])

    return line1, bar,


# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=np.linspace(
    1, 2*np.pi, 200), interval=1, blit=False)
plt.show()
