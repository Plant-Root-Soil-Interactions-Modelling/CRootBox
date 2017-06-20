import py_rootbox as rb
from rb_tools import *
import numpy as np
import matplotlib.pyplot as plt

rs = rb.RootSystem()
name = "Brassica_oleracea_Vansteenkiste_2014" 
rs.openFile(name)
rs.initialize() 

simtime = 60. # days
dt = 1.
N = round(simtime/dt) # steps
 
# Plot some scalar value over time
stype = rb.ScalarType.length 
stype_str = "length cm"
v_ = np.zeros(N)
v1_ = np.zeros(N)
v2_ = np.zeros(N)
v3_ = np.zeros(N)
for i in range(0,N):
    rs.simulate(dt,True)
    t = v2a(rs.getScalar(rb.ScalarType.type))
    v = v2a(rs.getScalar(stype))
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
plt.ylabel(stype_str)
plt.legend(["total", "tap root", "lateral", "2. order lateral"])
plt.savefig("results/example_3a.png")
plt.show()

# Find root tips and bases (two approaches)
rs.initialize() 
rs.simulate(4) # 4 days young....

print(rs.getNumberOfNodes(), "nodes")
print(rs.getNumberOfSegments(), "segments")

# Use polyline representation of the roots
polylines = rs.getPolylines()    
for r in polylines:
    print("base: \t", r[0])
    print("tip: \t", r[-1])  

# Or, use node indices to find tip or base nodes
nodes = vv2a(rs.getNodes())
tipI = rs.getRootTips()
baseI = rs.getRootBases()
for i in range(0,len(tipI)):
    print("base: \t", nodes[baseI[i]])
    print("tip: \t", nodes[tipI[i]])

