#!/usr/bin/python
# Simulation of 3x3 wheat root systems development in 120 days
# Distance betweem root systems is 40 cm
# Simulation output Root leng density [cm/cm3] (variable name RLD) and growing rate [cm/cm3/day] (variable name RLD_change) at each time step 
# created by Mai t.mai@jz-juelich.de
# 31 July 2017

import py_rootbox as rb
import numpy as np
import matplotlib.pyplot as plt

name = "wheat"
N = 3 # number of columns and rows
dist = 40 # distance between the root systems [cm]
Area = (dist*(N-1))**2 #[cm2]

# Creates and initializes N*N root systems
allRS = []
for i in range(0,N):
    for j in range(0,N):
         rs = rb.RootSystem()
         rs.openFile(name) 
         rs.getRootSystemParameter().seedPos = rb.Vector3d(dist*i,dist*j,-3) 
         allRS.append(rs)
         rs.initialize()

simtime = 120. # simulation time [day]
dt=1. # time step [day]
ntimestep= int(round(simtime/dt))
depth=150.; # soil depth [cm]
nl = 15; # number of layers

RLD_lastTime=vRLD=np.zeros((nl)) # Root length density at last time step [cm/cm3]
#simulation for each time step 
for t in range(0,ntimestep):
	for rs in allRS:
	    rs.simulate(dt)
	# Compute vertical RLD distribution in layers
	vRLD=np.zeros((N*N,nl)); # N*N rows, nl columns
	c=0
	for rs in allRS:
	      analysis = rb.SegmentAnalyser(rs)
	      rootLength = analysis.distribution(rb.ScalarType.length,0,depth,nl,True)
	      vRLD[c,:]=rootLength
	      vRLD[c,:] /= (depth/nl)# rootLength per soil depth
	      c += 1 # root system number
	RLD=np.sum(vRLD,axis=0)/(Area*(depth/nl))
	growingRate = (RLD-RLD_lastTime)/dt
	RLD_lastTime = RLD
	#Plot Root length density per layer
	z=np.linspace(0,depth*(-1),nl)   # depth*-1 is the (negativ) z coordinate
	plt.plot(RLD,z,'k-', color="blue",linewidth=2)
	x1,x2,y1,y2 = plt.axis()
	plt.axis((0,x2,y1,y2)) # set min of x axis to 0, because of some negative values of (mean-std)
	plt.xlabel('RLD (cm/cm3)')
	plt.ylabel('Depth (cm)')
	nameOutputFile="results/Soil3_wheat_RLD_time_"+str(t*dt)
	plt.savefig(nameOutputFile+".png")
	plt.close()
	#Plot Growing rate per layer
	plt.plot(growingRate,z,'k-', color="blue",linewidth=2)
	x1,x2,y1,y2 = plt.axis()
	plt.axis((0,x2,y1,y2)) # set min of x axis to 0, because of some negative values of (mean-std)
	plt.xlabel('Growing rate (cm/cm3/day)')
	plt.ylabel('Depth (cm)')
	nameOutputFile="results/Soil3_wheat_GrowingRate_time_"+str(t*dt)
	plt.savefig(nameOutputFile+".png")
	plt.close()

	# save data: z, mean, std in npz format (for multi array)
	#https://docs.scipy.org/doc/numpy/reference/generated/numpy.savez.html
	#np.savez(nameOutputFile, z=z, RLD=RLD)

# Export root systems (for visualisation) as single vtp files (as polylines)
c = 0
ana = rb.SegmentAnalyser() # see example 3b
for rs in allRS:
	c += 1 # root system number
	vtpname = "results/Soil3_wheat_"+str(c)+".vtp"
	rs.write(vtpname)
	ana.addSegments(rs) # collect all
	       
# Write all into single file (segments)
ana.write("results/Soil3_wheat_all.vtp") 

       
