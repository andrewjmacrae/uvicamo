import sys
from time import sleep
from daqhats import hat_list, HatIDs, mcc118
import numpy as np
from matplotlib import pyplot as plt

# get hat list of MCC daqhat boards
board_list = hat_list(filter_by_id = HatIDs.ANY)
if not board_list:
    print("No boards found")
    sys.exit()


board_num = 0
channel = 0

vlts = np.array([])
t = np.array([])
tTot = .1
dt = 0.001
nbr_smpls = int(tTot/dt)
tcurr = 0
board = mcc118(board_list[board_num].address)


for k in range(nbr_smpls):
    value = board.a_in_read(channel)    
    vlts = np.append(vlts,value)
    t = np.append(t,tcurr)
    tcurr += dt
    sleep(dt)

fig, ax = plt.subplots()

ax.plot(t,vlts,'-o')

plt.show()