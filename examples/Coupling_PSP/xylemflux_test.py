import math

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import time

import matplotlib.pyplot as plt

import py_rootbox as rb    
from rb_tools import *

from xylem_flux_ls import *

#
# initialize
#
rsname = "anagallis2010" 
rs = rb.RootSystem()
rs.openFile(rsname,"")
rs.initialize() # hydrotropism is not set right now, link to soil is missing

rs_Kr = ( 1e-10, 1e-10, 1e-10, 1e-10, 1e-10 )*10 # root hydraulic radial conductivity per root type 
rs_Kz = ( 1e-6, 1e-7, 1e-8, 1e-8, 1e-8) # root hydraulic axial conductivity per root type  

#
# simulate
#
rs.simulate(15)

#
# analysse
#
rs_nodes = rs.getNodes()
rs_seg = rs.getSegments()
rs_ana = rb.AnalysisSDF(rs) # segment analyser
rs_segL = vd2a(rs_ana.getScalar(rb.ScalarType.length))
rs_segType = vd2a(rs_ana.getScalar(rb.ScalarType.type))
rs_segR = vd2a(rs_ana.getScalar(rb.ScalarType.radius))

Q, b = xylem_flux_ls(rs_seg, rs_nodes, rs_segL, rs_segType, rs_segR, rs_Kr, rs_Kz)  
print(Q[0,1], Q[1,0])

s1 = [1,2] 
seg = [s1]

rho = 1e-3 # kg / cm      
g = 9.8 *100 # cm / s^2 
d = [-1000.*rho*g] 

Q, b = xylem_flux_bc_dirichlet(Q, b, seg, d)

x0 = d*np.ones(Q.shape[0]) # empirically proofen to be the best

t = time.time()
print("solve")
x, info = LA.cg(Q,b,x0=x0, tol=1e-12)  # tested with CG, CGS, GMRES, BICG. CG by far the best
print("fin")
print("CG: " + str(t-time.time()) +" sec" ) 

# plt.spy(Q, markersize=3) plt.show()
# print(Q)
# print(b)
print(x)
segX = nodes2seg(rs_nodes,rs_seg,x);
print(segX)
segX = a2vd(segX)
print(len(segX))

#
# output
# 
rs_ana.write(rsname+".vtp",segX)


