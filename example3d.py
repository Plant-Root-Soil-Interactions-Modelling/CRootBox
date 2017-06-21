import py_rootbox as rb
import math
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt

# sets all standard deviation to value*s
def set_all_sd(rs,s):
    for i in range(0,10):
        p = rs.getRootTypeParameter(i+1)
        p.lbs = p.lb*s
        p.las = p.la*s
        p.lns = p.ln*s
        p.nobs = p.nob*s
        p.rs = p.r*s
        p.a_s = p.a*s         

# Parameters
name = "Triticum_aestivum_a_Bingham_2011" 
simtime = 20
N = 100 # resolution of paramter
runs = 10 # iterations
theta0_ = np.linspace(0,math.pi/2,N)
dx_ = np.linspace(0.1,1.,N) # geometric 0.01 1?

# One simulation
def simulate(i):
    rs = rb.RootSystem()    
    rs.openFile(name)
    set_all_sd(rs,0.) # set all sd to zero                

    # vary parameter
    p1 = rs.getRootTypeParameter(1)
    p4 = rs.getRootTypeParameter(4)
    p1.theta = theta0_[i]
    p4.theta = theta0_[i]   

    # simulation
    rs.initialize()     
    rs.simulate(simtime, True) 

    # target 
    roots = rs.getPolylines()
    depth = 0.
    rad_dist = 0.
    for r in roots:
        depth += r[-1].z
        rad_dist += math.hypot(r[-1].x,r[-1].y)
        
    depth /= len(roots)
    rad_dist /= len(roots)

    return depth, rad_dist

# Parallel execution
param = [] # param is a list of tuples
for i in range(0,N):
    param.append((i,))
pool = Pool()
output =  pool.starmap(simulate,param)

depth_ = np.zeros(N)
rad_dist_ = np.zeros(N)
for i,o in enumerate(output):
    depth_[i] = o[0]
    rad_dist_[i] = o[1]
    
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10,8))
axes[0].set_xlabel('Theta (-)')
axes[1].set_xlabel('Theta (-)')
axes[0].set_ylabel('Mean tip depth (cm)')
axes[1].set_ylabel('Mean tip radial distance (cm)')
axes[0].plot(theta0_, depth_)
axes[1].plot(theta0_, rad_dist_)
fig.subplots_adjust()
plt.savefig("results/example_3d.png")
plt.show()

