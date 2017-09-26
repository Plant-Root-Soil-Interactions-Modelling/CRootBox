import py_rootbox as rb
from rb_tools import *
import math


def elongate(rs, inc, dt, se):

    ol = np.sum(v2a(rs.getScalar(rb.ScalarType.length)))
        
    rs_ = rb.RootSystem(rs) # copy
    rs_.simulate(dt, True)
    inc_ = np.sum(v2a(rs_.getScalar(rb.ScalarType.length))) - ol
    
    if inc_>inc and abs(inc_-inc)>0.5: # check if we have to perform a binary search  
        
        sl = 0. # left           
        sr = 1. # right
                        
        while abs(inc_-inc)>0.5: # binary search 
                                    
            m = (sl+sr)/2. # mid
            rs_ = rb.RootSystem(rs) # copy        
            se.setScale(m)  
            rs_.simulate(dt, True) 
            inc_ = np.sum(v2a(rs_.getScalar(rb.ScalarType.length))) - ol
            print("\tsl, mid, sr ", sl, m, sr, inc_)
            
            if inc_>inc: # concatenate
                sr = m
            else:
                sl = m                              
            
        return rs_                        
        
    else:
        return rs_
    
        

# Parameter
simtime = 30. # days
dt = 1
N = round(simtime/dt) # steps
maxinc = 20; # maximal length increment (cm/day), TODO base this value on some fancy model 

# Initialize root system
rs = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010" 
rs.openFile(name) 
rs.initialize() 

# Create proportional elongation callback 
se = rb.ProportionalElongation()
se.setScale(1.)
for i in range(0,10): # set scale elongation in the root type parameters
    p = rs.getRootTypeParameter(i+1)
    p.se = se # se = scale elongation 
 
ol = 0
 
# Simulation
for i in range(0,N):    
    print("\nSimulation ", i)
    rs = elongate(rs, maxinc*dt, dt, se)
    
    l = np.sum(v2a(rs.getScalar(rb.ScalarType.length)))
    inc =  l - ol
    ol = l
    print("elongated ", inc, " cm")
    
        
rs.write("results/example_carbon.vtp")



# 
# print("copy test")
# 
# rs = rb.RootSystem()
# name = "Anagallis_femina_Leitner_2010" 
# rs.openFile(name)     
# rs.initialize() 
# rs.simulate(20) # for a bit
# 
# rs2 = rb.RootSystem(rs) # copy the root system
# 
# nodes = vv2a(rs.getNodes())
# nodes2 = vv2a(rs2.getNodes())
# print(nodes.shape, nodes2.shape)
# 
# nodes2 = vv2a(rs2.getNodes())
# 
# rs.simulate(10)
# rs2.simulate(10)
# 
# nodes = vv2a(rs.getNodes())
# nodes2 = vv2a(rs2.getNodes())
# print(nodes.shape, nodes2.shape)
# 
# uneq = np.sum(nodes!=nodes2)
# print("Unequal nodes", uneq)


