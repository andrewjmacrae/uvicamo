import numpy as np
from Pi214 import *

setup()

pressed_stop = False

while pressed_stop == False:  # run until button is pressed

    if is_button_down(SW1):
        pressed_stop = True
        print("Button pressed - ending script")

    I = []
    theta = []
    
    tTotal = 1
    t = 0
    dt = 0.001

    #for now use an approximate theta
    freq = 10  # frequency of func gen
    theta = np.linspace(0, 2*np.pi*freq*tTotal, tTotal/dt)

    while t < tTotal:
        I.append(get_adc_val(scl=(5.0/1024)))

        sleep(dt)
        t+=dt
    
    #in the future we'll want Imax with no apparatus
    I = np.array(I)/np.max(I) #only want values between 0 and 1

    N = tTotal/dt

    A = (2/N)*np.sum(I)
    B = (4/N)*np.sum([I[i]*np.sin(2*theta[i]) for i in range(len(theta))])
    C = (4/N)*np.sum([I[i]*np.cos(4*theta[i]) for i in range(len(theta))])
    D = (4/N)*np.sum([I[i]*np.sin(4*theta[i]) for i in range(len(theta))])

    S0 = A-C
    S1 = 2*C
    S2 = 2*D
    S3 = B

    out = open('stokes_output.txt', 'w')
    out.write(str(S0))
    out.write(',')
    out.write(str(S1))
    out.write(',')
    out.write(str(S2))
    out.write(',')
    out.write(str(S3))
    out.close()

cleanup()