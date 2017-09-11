import py_rootbox as rb
from rb_tools import *
import math

# Parameter
simtime = 30. # days
dt = 1.
N = round(simtime/dt) # steps
maxinc = 20; # maximal length increment (cm/day), TODO base this value on some fancy model 

# Initialize root system
rs = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010" 
rs.openFile(name) 

# Set growth to linear (default is negative exponential)
for i in range(0,10):
    p = rs.getRootTypeParameter(i+1)
    p.gf = 2
     
rs.initialize() 

# Create proportional elongation callback 
se = rb.ProportionalElongation()
se.setScale(1.)

# Set scale elongation in the root type parameters
for i in range(0,10):
    p = rs.getRootTypeParameter(i+1)
    p.se = se # se = scale elongation 
 
ol = 0
se.setScale(1.) 

# Simulation
for i in range(0,N):
    
    se.setScale(1.) # rs and rs_ have the same se
     
    print("total length\t", ol)
    
    # calculate length increment
    rs_ = rb.RootSystem(rs) # copy
    rs_.simulate(dt, True)
    newl = np.sum(v2a(rs_.getScalar(rb.ScalarType.length)))
    inc = newl-ol

    print("unimpeded increment\t",inc)

    if inc>(maxinc*dt): # cm
        print("***")
        s = math.pow((maxinc*dt)/inc,0.6) # empirical, since root growth is not linear due to branching. smaller time steps could partially fix this
        se.setScale(s)    
             
    rs.simulate(dt, True) 
    
    l = np.sum(v2a(rs.getScalar(rb.ScalarType.length)))
    print("increment\t\t", l-ol)
    ol=l
    
    print()
    
rs.write("results/example_carbon.vtp")

