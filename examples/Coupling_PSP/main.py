import math

import numpy as np
from numpy import linalg as LA

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import PSP_infiltration1D as inf

import py_rootbox as rb    
from rb_tools import *

    
        
#
# Initialize soil domain
#    
isSuccess, soil = inf.readSoil("clay.txt")
if not isSuccess: 
    print("warning: wrong soil file.")
    quit()        
    
funcType = inf.CAMPBELL # inf.CAMPBELL, inf.RESTRICTED_VG, inf.IPPISCH_VG
        
Se = 0.5  # Initial degree of saturation ]0-1]

inf.initializeWater(funcType, soil, Se, inf.CELL_CENT_FIN_VOL)
        
# ubPotential = inf.airEntryPotential(funcType, soil[0]) # [J kg^-1] upper boundary condition INFILTRATION
ubPotential = inf.waterPotential(funcType, soil[0],inf.thetaFromSe(funcType, soil[0],0.1))    
isFreeDrainage = True

sumInfiltration = 0
totalIterationNr = 0

#
# Initialize root domain
#
rsname = "anagallis2010" 
rs = rb.RootSystem()
rs.openFile(rsname,"")
rs.initialize() # hydrotropism is not set right now, link to soil is missing
rs_Kr = ( 1e-10, 1e-10, 1e-10, 1e-10, 1e-10 ) # root hydraulic conductivity per root type 

#
# Simulation
#                   
simTime = 30*24 * 3600   
maxTimeStep = 24*3600                  
dt = 3600              
time = 0               

out_delay = 6*3600
out_next = out_delay             
out_img = np.zeros((inf.n+2,round(simTime/out_delay)))
out_c = 0

cflux = 0
  
#   
# Output
#
plt.ion()
fig = plt.figure() 

ax1 = fig.add_subplot(2, 2, 1)

ax2 = fig.add_subplot(2, 2, 2, projection='3d')  

ax3 = fig.add_subplot(2, 2, 3)
ax3.set_xlabel("Time [s]");
ax3.set_ylabel("Root uptake [1]") 

ax4 = fig.add_subplot(2, 2, 4)
ax4.set_ylim(-1, 0)
ax4.set_xlabel("Water content");
ax4.set_ylabel("Depth [m]") 



#
# Simulation loop
#
while (time < simTime):
    
    dt = min(dt, out_next - time) 

    #
    # Python Richards Code
    #    
    
    success, nrIterations, flux = inf.cellCentFiniteVolWater(funcType, soil, dt, ubPotential, isFreeDrainage,inf.LOGARITHMIC)
    #success, nrIterations, flux = (True,1,0)
    
    while (not success):
        print ("dt =", dt, "No convergence")
        dt = max(dt / 2, 1)
        for i in range(inf.n+2): # reset to old time step
            inf.theta[i] = inf.oldTheta[i]
            inf.psi[i] = inf.waterPotential(funcType, soil[inf.hor[i]], inf.theta[i])
        success, nrIterations, flux = inf.cellCentFiniteVolWater(funcType, soil, dt, ubPotential, isFreeDrainage,inf.LOGARITHMIC)
            
    for i in range(inf.n+2):
        inf.oldTheta[i] = inf.theta[i]
    sumInfiltration += flux * dt 
    time += dt
    
    #
    # Root System  
    #
    rs.simulate(dt/3600/24) # seconds to days
        
    #
    # Root System fluxes
    #    
    rs_nodes = rs.getNodes()
    rs_segI = rs.getSegments()
    rs_ana = rb.AnalysisSDF(rs) # segment analyser
    rs_segType = vd2a(rs_ana.getScalar(rb.ScalarType.type))
    rs_segL = vd2a(rs_ana.getScalar(rb.ScalarType.length))
    rs_segR = vd2a(rs_ana.getScalar(rb.ScalarType.radius))
    
    N = len(rs_segI)
    print("Number of segments is "+str(N))
    rs_segPot = np.ones(N) * (-1000) # pressure head
    
    rs_segV = np.zeros((N,3)) # direction of segment, needed for doussan
    rs_segZ = np.zeros(N) # segment mid point z coordinate
    for i in range(0,N):
        s = rs_segI[i]
        n1 = v2v(rs_nodes[s.x])
        n2 = v2v(rs_nodes[s.y])
        rs_segV[i,:] = n2-n1 # todo normalize
        mid = 0.5*(n1+n2)
        rs_segZ[i] = mid[2]  
     
    rs_flux = np.zeros(N)
    rs_flux2 = np.zeros(inf.n+2)
    for i in range(0,N):
        j = z2i(rs_segZ[i], inf.n)
        ind = int(rs_segType[i])
        rs_flux[i] =  (2*rs_segR[i]*math.pi*rs_segL[i])*rs_Kr[ind-1]*(inf.psi[j+1]-rs_segPot[i]) 
        #print((inf.psi[j+1]-rs_segPot[i]))
        rs_flux2[j+1] += max((rs_flux[i]*dt),0)             
    #print(rs_flux2)                  
    #
    # Apply Sink
    #    
    for i in range(0,inf.n):
        inf.theta[i+1] -= ( rs_flux2[i+1]/inf.area/inf.dz[i+1] )    
        inf.theta[i+1] = min(inf.theta[i+1],1)
        inf.theta[i+1] = max(inf.theta[i+1],0)          
      
    print(inf.theta)
    for i in range(0,inf.n+2):
        inf.psi[i] = inf.waterPotential(funcType, soil[inf.hor[i]], inf.theta[i])
     
    cflux += float(np.sum(rs_flux2))   
     
    #
    # Output
    #
    print("time =", int(time), "\tdt =", dt, "\tIter. =", int(nrIterations),  "\tInf:", format(sumInfiltration, '.3f'))
    
    if time>=out_next:
        out_next += out_delay
        
        # Subplot 1
        out_img[:,out_c] = inf.theta[:]
#         imin = np.min(np.min(out_img))
#         imax = np.max(np.max(out_img))     
        out_c += 1
        #ax1.imshow((out_img-imin)/(imax-imin), interpolation="nearest", cmap="hot")        
        ax1.imshow(out_img, interpolation="nearest", cmap="hot")        

        # Subplot2
        plotRSscatter(ax2, rs.getRootTips())        
        
        # Subplot 3                    
        ax3.plot([time,out_next], [cflux,cflux], 'k-')
    
        # Subplot 4
        ax4.plot(inf.theta,-inf.z,'k-')   
    
        plt.pause(0.0001)
            
    if (float(nrIterations/inf.maxNrIterations) < 0.1): 
        dt = min(dt*2, maxTimeStep) # go fasetr        
           
plt.ioff()
plt.show()
