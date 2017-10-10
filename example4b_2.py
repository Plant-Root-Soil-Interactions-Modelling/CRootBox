import py_rootbox as rb
import math
import numpy as np
from rb_tools import *

rs = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010" 
rs.openFile(name)

soilprop = rb.SoilProperty1Dlinear(0,-50,100) # for root elongation 
data = np.ones((100,))
soilprop.setData(a2v(data))
  
# print("Value at -2 ", soilprop.getValue(rb.Vector3d(0,0,-2), r))
  
# Manually set scaling function and tropism parameters
sigma = [0.4, 1., 1., 1., 1. ] * 2
for i in range(0,10):  
    p = rs.getRootTypeParameter(i+1)
    p.dx = 0.25 # adjust resolution
    p.tropismS = sigma[i] 
     
    # 1. Scale elongation
    p.se = soilprop
       
# simulation
rs.initialize()
simtime = 120.
dt = 1.
N = 120/dt
for i in range(0,round(N)):
    rs.simulate(dt,True)


print("fin")