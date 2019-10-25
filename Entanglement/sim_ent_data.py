import numpy as np
from time import sleep


dx = 0.05
fileName = 'spdcdata.txt'


def pull(fName):
    pullData = open(fName,'r').read()
    dataArray = pullData.split('\n')
    det1=[]
    det2=[]
    coinc=[]
    for eachLine in dataArray:
        if len(eachLine)>1:
            d1,d2,c2x = eachLine.split(',')
            det1.append(float(d1))
            det2.append(float(d2))
            coinc.append(float(c2x))
    return det1,det2,coinc

def push(fName,d1,d2,c2x,d10,d20,c2x0,nMax):

    if len(d1) < nMax-1:        
        d1.append(d10)
        d2.append(d20)
        c2x.append(c2x0)
    else:
        d1[0:-1] = d1[1:]
        d2[0:-1] = d2[1:]
        c2x[0:-1] = c2x[1:]
        
        d1[-1] = d10
        d2[-1] = d20
        c2x[-1] = c2x0

    with open(fName,'w+') as f:
        for k in range(len(d1)):
            f.write(str(d1[k])+','+str(d2[k])+','+str(c2x[k])+'\n')

def init(fName):
    with open(fName,'w+') as f:
        f.write('0,0,0\n')



init(fileName)

t=0.0
Ncnt = 1e3

def acc2x(N1,N2,dT,mm):
    tTot = 1.
    M = tTot/dT
    C1 = np.zeros(int(M))
    C2 = np.zeros(int(M))
    for k in range(int(N1)):
        C1[int(np.random.rand()*M)] = 1
    for k in range(int(N2)):
        C2[int(np.random.rand()*M)] = 1
    Coinc = C1*C2

    return sum(Coinc)+mm*min(N1,N2)

while True:
    y1 = Ncnt*(np.random.rand()*.25+1.5+((np.sin(.01*t)+np.sin(.04*t)+np.sin(.07*t)))/3)/2.5
    y2 = Ncnt*(np.random.rand()*.25+1.5+((np.cos(.02123*t)+np.sin(.05114*t)+np.sin(.061237*t)))/3)/2.5
    y3 = acc2x(y1,y2,2e-6,0.2)
    det1,det2,coinc = pull(fileName)
    
    push(fileName,det1,det2,coinc,y1,y2,y3,50)
    sleep(.25)
    t+=1.7