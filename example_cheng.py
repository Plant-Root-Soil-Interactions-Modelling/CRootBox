#
# Example Exudation
#
# 1) Opens plant and root parameters from a file
# 2) Simulates root growth
# 3) Outputs a VTP (for vizualisation in ParaView)
#
#  Computes analytical solution of moving point/line sources based on Carslaw and Jaeger
#

import py_rootbox as rb
import numpy as np
# from evtk.hl import gridToVTK # how to install? TVTK ?


def v2a(vd): # rb.std_vector_double_ to numpy array    
    l = np.zeros(len(vd)) 
    for i in range(0,len(vd)):
        l[i] = vd[i]
    return l





rootsystem = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010" 

#
# Open plant and root parameter from a file
#
rootsystem.openFile(name)

#
# Initialize
#
rootsystem.initialize() 

#
# Simulate
#
simtime = 30  # or 20, 40, 60 days
dt = 10 # try other values here
N = round(simtime/dt) # steps
for i in range(0,int(N)):
    rootsystem.simulate(dt);

#
# Export final result (as vtp)
#
rootsystem.write("results/"+name+".vtp",rb.OutputType.segments) # use ot_polylines for nicer visualization, ot_segments for animations


params = rb.ExudationParameters()  

X = 10
Y = 10
Z = 100
width = 30
depth = 30

C = rb.getExudateConcentration(rootsystem, params, 10, 10, 100, 30, 50);
C = v2a(C); # make a numpy array 
C = np.reshape(C, (10,10,100) ) # hope that works




