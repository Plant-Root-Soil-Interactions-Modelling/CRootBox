#
# The Python version of example1.h:
#
#
# 1) Opens plant and root parameters from a file
# 2) Simulates root growth
# 3) Outputs a VTP (for vizualisation in ParaView)
#    In Paraview: use tubePLot.py script for fancy visualisation (Macro/Add new macro...), apply after opening file
#
#  Additionally, exports the line segments as .txt file to import into Matlab for postprocessing
#
import py_rootbox as rb

rootsystem = rb.RootSystem()
name = "anagallis2010" 

#
# Open plant and root parameter from a file
#
rootsystem.openFile(name,"modelparameter/")
# rootsystem.writeParameters() # not exposed to python yet

#
# Initialize
#
rootsystem.initialize()

#
# Simulate
#
simtime = 30  # or 20, 40, 60 days
dt = 1 # try other values here
N = round(simtime/dt) # steps

for i in range(0,int(N)):
    rootsystem.simulate(dt);

#
# Export final result (as vtp)
#
rootsystem.write(name+".vtp",rb.OutputType.segments) # use ot_polylines for nicer visualization, ot_segments for animations

#
# Export segments for Matlab analysis
#    
analysis = rb.AnalysisSDF(rootsystem)
analysis.write(name+".txt")

#
# Total length and surface
#     
l = analysis.getSummed(rb.ScalarType.length)  
print("Total root system length is "+str(l)+" cm")    

print("Finished with a total of "+str(rootsystem.getNumberOfNodes())+ " nodes")

# end of example 1  

#
# Python testing ... 
#
L = analysis.getScalar(rb.ScalarType.length); # we can get scalar Data and do soemthing with it :-)
lt= 0
for l in L :
    lt += l
print(lt)


rtp = rb.RootTypeParameter()
p = rtp.realize();
print(p)


ln_ = rb.std_vector_double_(); # how to append/fill values????
#for i in range(0,3):
    #ln_[i]=0.1;
    
from math import inf
    
type =1 # 1
lb = 1 # cm
la = 7  # cm
nob = 0 # 1
r = 2 #  cm/day
a = 0.1 # cm
theta = 1.2 # rad
rlt = inf # day

p2 = rb.RootParameter()
p2.set(type,lb, la,ln_, nob, r, a, theta, rlt)
print(p2)







