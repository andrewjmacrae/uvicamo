import numpy as np

I_bg = 56e-3
I_min = .92
I_max = 3.60

eta = (I_min-I_bg)/(I_max-I_bg)

sig = 2*eta-1
print(f'retardance = arccos({round(sig,3)}) = {round(np.arccos(sig),3)}')
