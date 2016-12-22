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
import math
import random
import numpy as np
import matplotlib.pyplot as plt


rootsystem = rb.RootSystem()
name = "anagallis2010"
nameOutputFile="Mirjam"

allRS = [ ]

#
# Creates N root systems
#
N=100
for i in range(0,N-1):
     rs = rb.RootSystem()
     rs.openFile(name) 
     allRS.append(rs)

#
# Simulate
#
simtime = 28
for rs in allRS:
    rs.setSeed(random.random() )
    rs.initialize()
    rs.simulate(simtime)
    s=rs.getSegments()
    print(s)

#
# Export results as single vtp files
#
c = 0
for rs in allRS:
      c += 1 # root system number
      vtpname = "results/"+name+str(c)+".vtp";
      rs.write(vtpname, rb.OutputType.polylines);

#
# Compute vertical RLD distribution in layers 
#   
nl = 100; # number of layers
vRLD=np.zeros((N,nl)); # N rows, nl columns
depth=100;
c=0
for rs in allRS:
      c += 1 # root system number
      analysis = rb.AnalysisSDF(rs)
      RLD = analysis.distribution(rb.ScalarType.length,0,depth,nl,1)
      vRLD[c,:]=RLD

z=np.linspace(0,depth*(-1),nl)   # depth*-1 is the (negativ) z coordinate 
mean=np.mean(vRLD,axis=0)
std=np.std(vRLD,axis=0)
plt.figure(figsize=(3.8,3))
plt.plot(mean,z,'k-', color="blue",linewidth=2)
plt.fill_betweenx(z,mean+std,mean-std,color="blue",edgecolor='',alpha=0.5)
x1,x2,y1,y2 = plt.axis()
plt.axis((0,x2,y1,y2)) # set min of x axis to 0, because of some negative values of (mean-std)
plt.savefig(nameOutputFile+".png")
plt.show()
# save data: z, mean, std in npz format (for multi array)
#https://docs.scipy.org/doc/numpy/reference/generated/numpy.savez.html
np.savez(nameOutputFile, z=z, mean=mean, std=std)
