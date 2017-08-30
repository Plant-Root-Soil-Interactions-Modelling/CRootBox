#PSP_plant.py
from __future__ import division
from PSP_public import *
from PSP_soil import *

t = np.zeros(n+2)
Rr = np.zeros(n+2)
rs = np.zeros(n+2)
bz = np.zeros(n+2)
leafPot=np.zeros(n+2)
#pb=np.zeros(n+2)
    
def InitTransp(soil, z):
    rw = 25000000000
    r1 = 0.001
    for i in range(1, n+1):
        if (z[i] > rootMin and z[i] < rootDepth):
            L = 40000 * (rootDepth - z[i]) / rootDepth
            Rr[i] = 2 * rw / (L * (z[i+1] - z[i-1]))
            bz[i] = ((1 - soil.Campbell_n) * 
            np.log(np.pi * r1 * r1 * L) / (2 * np.pi * L * (z[i+1] - z[i-1])))
        else:
            Rr[i] = 1E+20
            bz[i] = 0
    return

def PlantWaterUptake(tp, funcType, soil, theta, psi, z):
    pc = -1500
    RL = 3000000
    sp = 10
    pl = 0
    pb = 0
    rb = 0
    for i in range(1, n+1):
        rs[i] = bz[i] / hydraulicConductivityFromTheta(funcType, soil, theta[i])
        pb = pb + psi[i] / (Rr[i] + rs[i])
        rb = rb + 1 / (Rr[i] + rs[i])
    pb = pb / rb
    rb = 1 / rb
    if (pl > pb):
        pl = pb - tp * (RL + rb)
    f = 11
    while(abs(f)>10):
        xp = np.power(pl / pc,sp)
        sl = tp * (RL + rb) * xp * sp / (pl * (1 + xp) * (1 + xp)) - 1
        f = pb - pl - tp * (RL + rb) / (1 + xp)
        pl = pl - f / sl  
    tr = tp / (1 + xp)
    for i in range(1, n+1):  
        t[i] = (psi[i] - pl - RL * tr) / (Rr[i] + rs[i])  
    return tr,pl,pb
