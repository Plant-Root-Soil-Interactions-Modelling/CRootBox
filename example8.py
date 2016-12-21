#
# Example 7
#
# 100 root system
# Chemotropism ,
# proof of concept, with a static nutrient concentration
#
import py_rootbox as rb
import math
import random
import numpy as np
import matplotlib.pyplot as plt

name = "maize"

#
# Plant and root parameter from a file
#
allRS = [ ]
allRSnormal = [ ]
#
# Creates N root systems
#
N=100
for i in range(0,N-1):
	rootsystem1 = rb.RootSystem();
	rootsystem1.openFile(name);
	allRS.append(rootsystem1);
	rootsystem2 = rb.RootSystem();
	rootsystem2.openFile(name);
	allRSnormal.append(rootsystem2)

# Created soil property
#
sideBox = rb.SDF_PlantBox(30,30,100)
#rootsystem.setGeometry(sideBox)  # for vizualisation

#nutrientBox = rb.SDF_PlantBox(30,30,5)
#nutrientBox2 = rb.SDF_RotateTranslate(nutrientBox, rb.Vector3d(0,0,-20))
nutrientBox = rb.SDF_PlantBox(10,10,10)
nutrientBox2 = rb.SDF_RotateTranslate(nutrientBox, rb.Vector3d(-5,0,-10))


maxS = 0.7 # maximal saturation
minS = 0.1 # minimal saturation
slope = 20 # [cm] linear gradient between min and ma
soilprop = rb.SoilPropertySDF(nutrientBox2, maxS, minS, slope)

# set the soil properties before calling initialize
#rootsystem.setSoil(soilprop);

#
# Initialize
#
#rootsystem.initialize(4,5); # it is important to call initialize() after setGeometry()

#
# Simulate
#
simtime = 60; # e.g. 30 or 60 days
for rs in allRS:
	rs.setGeometry(sideBox)
	for i in range(0,10):
	    rs.getRootTypeParameter(i+1).tropismT = rb.TropismType.chemo;
	    rs.getRootTypeParameter(i+1).tropismN = 2; # N
	    rs.getRootTypeParameter(i+1).tropismS = 0.5; # sigma
	print (rs.getRootSystemParameter().seedPos)
	rs.setSoil(soilprop)
	rs.setSeed(random.random())
	rs.initialize()
	rs.simulate(simtime)
	s=rs.getSegments()
	print(s)

for rs in allRSnormal:
	rs.setGeometry(sideBox)
	rs.setSeed(random.random())
	rs.initialize()
	rs.simulate(simtime)
	s_normal=rs.getSegments()
	print(s_normal)
#
# Export results as single vtp files
#
c = 0
for rs in allRS:
      c += 1 # root system number
      vtpname = "results/"+name+str(c)+".vtp";
      rs.write(vtpname, rb.OutputType.polylines);
c = 0
for rs in allRSnormal:
      c += 1 # root system number
      vtpname = "results/"+name+"_normal"+str(c)+".vtp";
      rs.write(vtpname, rb.OutputType.polylines);

#
# Export container geometry as Paraview Python script (run file in Paraview by Tools->Python Shell, Run Script)
#
rs.write(name + "_example8.py",0);

#
# Compute vertical RLD distribution in layers
#
nl = 50; # number of layers
vRLD=np.zeros((N,nl)); # N rows, nl columns
depth=100;
c=0
for rs in allRS:
      c += 1 # root system number
      analysis = rb.AnalysisSDF(rs)
      RLD = analysis.distribution(rb.ScalarType.length,0,depth,nl,1)
      vRLD[c,:]=RLD

vRLD_normal=np.zeros((N,nl)); # N rows, nl columns
c=0
for rs in allRSnormal:
      c += 1 # root system number
      analysis = rb.AnalysisSDF(rs)
      RLD = analysis.distribution(rb.ScalarType.length,0,depth,nl,1)
      vRLD_normal[c,:]=RLD

z=np.linspace(0,depth*(-1),nl)   # depth*-1 is the (negativ) z coordinate
mean=np.mean(vRLD,axis=0)
std=np.std(vRLD,axis=0)
mean_normal=np.mean(vRLD_normal,axis=0)
std_normal=np.std(vRLD_normal,axis=0)
plt.figure(figsize=(3.8,3))
plt.plot(z,mean,'k-',linewidth=2)
plt.fill_between(z,mean+std,mean-std,color='#b9cfe7',edgecolor='')
plt.plot(z,mean_normal,'ko',color="red")
plt.fill_between(z,mean_normal+std_normal,mean_normal-std_normal,color='red',edgecolor='')
plt.show()

#print("Finished with a total of " + str(rootsystem.getNumberOfNodes()) + " nodes\n")
