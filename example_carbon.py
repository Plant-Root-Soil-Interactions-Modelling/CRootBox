import py_rootbox as rb
from rb_tools import *

# Initialize root system
rs = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010" 
rs.openFile(name) 
rs.initialize() 

# Parameter
simtime = 60. # days
dt = 1.
N = round(simtime/dt) # steps

# Simulation
for i in range(0,N):
    
    v = v2a(rs.getScalar(rb.ScalarType.volume))
    vol = np.sum(v)
    
    # 
    rs_ = rb.RootSystem(rs) # copy
    rs_.simulate(dt)
    v = v2a(rs_.getScalar(rb.ScalarType.volume))
    print(len(v))
    newvol = np.sum(v)
    print(vol, newvol)

    #
    
    vinc = newvol-vol;
    # print("\n", vinc)
    
    rs.simulate(dt, True) 






