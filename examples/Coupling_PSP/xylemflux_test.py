import math
import time
import os # add the search path for py_rootbox.so (probably there is a nicer way to do it?)
import sys
cwd = os.getcwd()
i = cwd.index("CRootBox")
sys.path.append(cwd[0:i+8]) 

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA

import matplotlib.pyplot as plt

import py_rootbox as rb    
from rb_tools import *

import xylem_flux 

#
# initialize (root system)
#
path = parameterPath()
rsname = "anagallis_Leitner_et_al(2010)" 
rs = rb.RootSystem()
rs.openFile(rsname,path)
rs.initialize() # hydrotropism is not set right now, link to soil is missing

#
# simulate
#
rs.simulate(15)

#
# results 
#
seg = seg2a(rs.getSegments())
nodes = vv2a(rs.getNodes())
rs_ana = rb.SegmentAnalyser(rs) 
type = v2a(rs_ana.getScalar(rb.ScalarType.type))
radius = v2a(rs_ana.getScalar(rb.ScalarType.radius))

#
# initialize (xylem_flux)
#
rs_Kr = np.array([ 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10 ]) # root hydraulic radial conductivity per root type [cm^2 s / kg], plausible values?
rs_Kz = np.array([ 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8 ]) # root hydraulic axial conductivity per root type [cm^5 s], plausible values?  
kr = np.array(list(map(lambda t: rs_Kr[int(t)-1], type))) # convert from 'per type' to 'per segment'
kz = np.array(list(map(lambda t: rs_Kz[int(t)-1], type)))

rho = 1e-3 # kg / cm      
g = 9.8 *100 # cm / s^2 

soil_p = lambda x,y,z : -60 # kg/(cm s^2)

#
# create linear system
#
Q, b = xylem_flux.linear_system(seg, nodes, radius, kr, kz, rho, g, soil_p)  
# print("row 0 ")
# print(Q[0,0], Q[0,1])
# print(Q[1,0], Q[1,1])
# print(b[0])
# print(b[1])

#
# apply Dirichlet BC
#
d = [-1000.*rho*g] 
n0= np.array([0]) # node indices
# Q, b = xylem_flux.bc_dirichlet(Q, b, n0, d)

#
# apply Neumann BC
#
seg0 = np.array([0]) # segment indices
aflux = np.array([-2.e-6])
Q, b = xylem_flux.bc_neumann(Q, b, seg0, aflux, seg, nodes)

#
# solve LS
#
t = time.time()
print("solve")
# x0 = d*np.ones(Q.shape[0]) # empirically proofen to be the best
# x, info = LA.cg(Q,b,x0=x0, tol=1e-12)  # tested with CG, CGS, GMRES, BICG. CG by far the best
x = LA.spsolve(Q, b) # direct
print("fin")
# print("CG: " + str(time.time()-t) +" sec" ) 
print("spsolve: " + str(time.time()-t) +" sec" )

#
# output
# 
segP = nodes2seg(nodes,seg,x) 
axial_flux = xylem_flux.axial_flux(x, seg, nodes, kz, rho, g)
radial_flux = xylem_flux.radial_flux(x, seg, nodes, radius, kr, soil_p)
net_flux = axial_flux+radial_flux

rs_ana.addUserData(a2v(segP),"pressure")
rs_ana.addUserData(a2v(axial_flux),"axial_flux")
rs_ana.addUserData(a2v(radial_flux),"radial_flux")
rs_ana.addUserData(a2v(net_flux),"net_flux")
rs_ana.write(rsname+".vtp")




