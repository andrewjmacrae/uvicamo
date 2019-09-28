import numpy as np 

#get I 

#get theta 

A = (2/N)*np.sum(I)
B = (4/N)*np.sum([I[i]*np.sin(2*theta[i]) for i in range(len(theta))])
C = (4/N)*np.sum([I[i]*np.cos(4*theta[i]) for i in range(len(theta))])
D = (4/N)*np.sum([I[i]*np.sin(4*theta[i]) for i in range(len(theta))])

S0 = real(A-C)
S1 = real(2*C)
S2 = real(2*D)
S3 = real(B)

S = [S0,S1,S2,S3]


#out = open('stokes_output.txt','w')
#out.write(S)
