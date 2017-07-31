import scipy.sparse.linalg as LA
from scipy import sparse

import py_rootbox as rb    
from rb_tools import *

import xylem_flux 

# Simulate a root system
name = "Sorghum_bicolor_NA_NA"
simtime = 60 # days
dt = 0.1 
N = round(simtime/dt)
rs = rb.RootSystem()
rs.openFile(name)
rs.initialize() 

# Incrementally build an adjacency matrix 
nodes = np.array([])
seg= np.array([])


print()
for i in range(0,N):    
     
     rs.simulate(dt) 
     
     print("Number of nodes", rs.getNumberOfNodes()) # equals the number of new segments
     print("Number of roots", len(rs.getRootTips()))
     
     print("Number of new nodes", rs.getNumberOfNewNodes()) # equals the number of new segments
     print("Number of new roots", rs.getNumberOfNewRoots())

     uni = v2a(rs.getUpdatedNodeIndices())    
     unodes =  vv2a(rs.getUpdatedNodes())
     print("Number of node updates", len(unodes), len(uni))
             
     newnodes = vv2a(rs.getNewNodes())
     newsegs = seg2a(rs.getNewSegments())
     
     print()

