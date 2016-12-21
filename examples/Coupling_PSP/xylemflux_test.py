import math

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import py_rootbox as rb    
from rb_tools import *
from xylem_flux import *

#
# initialize
#
rsname = "anagallis2010" 
rs = rb.RootSystem()
rs.openFile(rsname,"")
rs.initialize() # hydrotropism is not set right now, link to soil is missing

rs_Kr = ( 1e-10, 1e-10, 1e-10, 1e-10, 1e-10 ) # root hydraulic radial conductivity per root type 
rs_Kz = ( 1e-8, 1e-8, 1e-8, 1e-8, 1e-8 ) # root hydraulic axial conductivity per root type  

#
# simulate
#
rs.simulate(1)

#
# analysse
#
rs_nodes = rs.getNodes()
rs_segI = rs.getSegments()
rs_ana = rb.AnalysisSDF(rs) # segment analyser
rs_segL = vd2a(rs_ana.getScalar(rb.ScalarType.length))
rs_segType = vd2a(rs_ana.getScalar(rb.ScalarType.type))
rs_segR = vd2a(rs_ana.getScalar(rb.ScalarType.radius))

Q, b = xylem_flux(rs_segI, rs_nodes, rs_segL, rs_segType, rs_segR, 0)

print("solve")
x, info = LA.bicg(Q,b,x0=None, tol=1e-16)
print("fin")

plt.spy(Q, markersize=3)
plt.show()
print(Q)
print(b)
print(x)


#
# output
# 
rs_ana.write(rsname+".vtp")


