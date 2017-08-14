import py_rootbox as rb
from rb_tools import *
import math
import numpy as np

# Parameter
simtime = 30. # days
dt = 1.
N = round(simtime/dt) # steps
name = "Anagallis_femina_Leitner_2010" 

print("\n1 original")
rs = rb.RootSystem()
rs.openFile(name) 
rs.setSeed(1.) # before initialize to mimic
rs.initialize()
rs2 = rb.RootSystem(rs)
rs.simulate(simtime)
print("total length", np.sum(v2a(rs.getScalar(rb.ScalarType.length))))

print("\n2 copy")
rs2.simulate(simtime)
print("total length", np.sum(v2a(rs2.getScalar(rb.ScalarType.length))))

print("\n3 rebuild same")
rs3 = rb.RootSystem()
rs3.openFile(name) 
rs3.setSeed(1.)
rs3.initialize()
rs3.simulate(simtime)
print("total length", np.sum(v2a(rs3.getScalar(rb.ScalarType.length))))

    
    
    
       