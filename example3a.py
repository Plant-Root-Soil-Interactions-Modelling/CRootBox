import py_rootbox as rb
from rb_tools import *

import numpy as np
import matplotlib.pyplot as plt

rootsystem = rb.RootSystem()
name = "Brassica_oleracea_Vansteenkiste_2014" # Brassica_oleracea_Vansteenkiste_2014, Brassica_napus_a_Leitner_2010
rootsystem.openFile(name)
rootsystem.initialize() 

simtime = 60. # days
dt = 1.
N = round(simtime/dt) # steps
 
#
# plot some scalar value over time
#
type = rb.ScalarType.length #  type, radius, order, time, length, surface, one, parenttype
type_str = "length cm"
v_ = np.zeros(N)
v1_ = np.zeros(N)
v2_ = np.zeros(N)
v3_ = np.zeros(N)
for i in range(0,N):
    rootsystem.simulate(dt,True)
    t = v2a(rootsystem.getScalar(rb.ScalarType.type))
    v = v2a(rootsystem.getScalar(type))
    v_[i] = np.sum(v)
    v1_[i] = np.sum(v[t==1])  
    v2_[i] = np.sum(v[t==2]) 
    v3_[i] = np.sum(v[t==3]) 
        
t_ = np.linspace(dt, N*dt, N)
plt.plot(t_,v_)
plt.plot(t_,v1_)
plt.plot(t_,v2_)
plt.plot(t_,v3_)
plt.xlabel("time (days)")
plt.ylabel(type_str)
plt.legend(["total", "tap root", "lateral", "2. order lateral"])
plt.savefig("results/example_3a.png")
plt.show()

#
# find root tips and bases (two approaches)
#
rootsystem.initialize() 
rootsystem.simulate(4) # 4 days young....

print(rootsystem.getNumberOfNodes(), "nodes")
print(rootsystem.getNumberOfSegments(), "segments")

# use polyline representation of the roots
polylines = rootsystem.getPolylines()    
for r in polylines:
    print("base: \t", r[0])
    print("tip: \t", r[-1])  

# or use node indices to find tip or base nodes
nodes = vv2a(rootsystem.getNodes())
tipI = rootsystem.getRootTips()
baseI = rootsystem.getRootBases()
for i in range(0,len(tipI)):
    print("base: \t", nodes[baseI[i]])
    print("tip: \t", nodes[tipI[i]])

