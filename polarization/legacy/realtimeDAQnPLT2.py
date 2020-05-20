import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#define functions
def polarization_ellipse(S):
    '''
    given a stokes vector, this function plots the corresponding
    polarization ellipse

    returns a list of x and y values on the Ex-Ey plane
    '''

    #set individual stokes parameters from stokes vector
    S0, S1, S2 = S[0], S[1], S[2]

    #solve for psi, the angle from the x axis of the ellipse
    if S1 == 0:
        psi = np.pi/4
    else:
        psi = 0.5*np.arctan(S2/S1)

    #define ellipse parameters from stokes vectors
    a = np.sqrt(0.5*(S0+np.sqrt(S1**2+S2**2)))
    b = np.sqrt(0.5*(S0-np.sqrt(S1**2+S2**2)))
    rot = np.matrix([[np.cos(psi), -1*np.sin(psi)],
                     [np.sin(psi), np.cos(psi)]])
    ba = b/a
    x1, x2, y1, y2 = [], [], [], []
    #create an x array for plotting ellipse y values
    x = np.linspace(-a, a, 200)

    for x in x:
        #cartesian equation of an ellipse
        Y1 = ba*np.sqrt(a**2-x**2)
        #Y1 relection about the x-axis
        #rotate the ellipse by psi
        XY1 = np.matrix([[x],
                         [Y1]])
        XY2 = np.matrix([[x],
                         [-Y1]])
        y1.append(float((rot*XY1)[1]))
        x1.append(float((rot*XY1)[0]))
        y2.append(float((rot*XY2)[1]))
        x2.append(float((rot*XY2)[0]))

    #x2,y2 reversed in order so that there is continuity in the ellipse (no line through the middle)
    x = x1+x2[::-1]
    y = y1+y2[::-1]

    return x, y

def animate(S):

    # animation function.  This is called sequentially
    #later, iterate through S instead of phi

    x1, y1 = [], []

    print('in animate, S = ', S)

    x1, y1 = polarization_ellipse(S)

    line1.set_data(x1, y1)
    for i in range(len(S)):
        bar[i].set_height(S[i])

    return line1, bar,

def get_stokes():
    N = 1
    I = []
    theta = []
    while True:
        #collect data from adc
        #I = hat.a_in_scan_read(read_request_size, timeout)
        #N = read_request_size
        #theta = fake_theta

        #fake data for now
        N += 1.0
        theta_tmp = N/5.0*np.pi
        I.append(np.cos(theta_tmp)**2)
        #I.append((np.cos(theta_tmp)+np.cos(4*theta_tmp)+1.9)/3.5)
        #I.append(0.9)
        theta.append(theta_tmp)

        A = (2.0/N)*np.sum(I)
        B = (4.0/N)*np.sum([I[i]*np.sin(2*theta[i])
                            for i in range(len(theta))])
        C = (4.0/N)*np.sum([I[i]*np.cos(4*theta[i])
                            for i in range(len(theta))])
        D = (4.0/N)*np.sum([I[i]*np.sin(4*theta[i])
                            for i in range(len(theta))])

        #returns [S0,S1,S2,S3]
        yield [np.real(A-C), np.real(2.0*C), np.real(2.0*D), np.real(B)]

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

#initialize bar graph as a variable outside animation function so
#  we don't have to clear out data each frame
bar = plt.bar([0, 1, 2, 3], [0, 0, 0, 0], align='center')

anim = animation.FuncAnimation(fig, animate, get_stokes, interval=100)
plt.show()