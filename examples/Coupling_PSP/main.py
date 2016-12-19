#PSP_infiltrationRedistribution1D
from __future__ import print_function, division

import matplotlib.pyplot as plt
import PSP_infiltration1D as inf

import py_rootbox as rb    
import numpy as np


#
#    auxiliary functions that should be moved to py_rootbox
#
def v2v(v): # rb.Vector3 to numpy array 
    return np.array([v.x, v.y, v.z])

def vd2l(vd): # rb.std_vector_double_ to numpy array
    N  = len(vd)
    l = np.zeros(N) 
    for i in range(0,N):
        l[i] = vd[i]
    return l
        
#
# Initialize soil domain
#    
isSuccess, soil = inf.readSoil("soilUniform.txt")
if not isSuccess: 
    print("warning: wrong soil file.")
    quit()        
    
funcType = inf.CAMPBELL # inf.CAMPBELL, inf.RESTRICTED_VG, inf.IPPISCH_VG
        
Se = 0.5  # Initial degree of saturation ]0-1]

inf.initializeWater(funcType, soil, Se, inf.CELL_CENT_FIN_VOL)
        
ubPotential = inf.airEntryPotential(funcType, soil[0]) # [J kg^-1] upper boundary condition
#ubPotential = inf.waterPotential(funcType, soil[0],inf.thetaFromSe(funcType, soil[0],0.1))    
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

#
# Simulation
#                   
simTime = 31 * 24 * 3600   
maxTimeStep = 24*3600                  
dt = 3600              
time = 0                            
  
#   
# Output
#
plt.ion()
f, myPlot = plt.subplots(2, figsize=(10, 8), dpi=80)
myPlot[1].set_xlim(0, simTime)
myPlot[1].set_ylim(0, 0.1)
myPlot[1].set_xlabel("Time [s]",fontsize=16,labelpad=8)
myPlot[1].set_ylabel("Infiltration Rate [kg m$^{-2}$ s$^{-1}$]",fontsize=16,labelpad=8)
plt.tick_params(axis='both', which='major', labelsize=12,pad=6)

#
# Simulation loop
#
while (time < simTime):
    dt = min(dt, simTime - time) 

    #
    # Python Richards Code
    #    
    success, nrIterations, flux = inf.cellCentFiniteVolWater(funcType, soil, dt, ubPotential, isFreeDrainage,inf.LOGARITHMIC)
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
    rs_segType = vd2l(rs.getScalar(rb.OutputType.segments, rb.ScalarType.type))
    rs_segL = vd2l(rs.getScalar(rb.OutputType.segments, rb.ScalarType.length))
    rs_segR = vd2l(rs.getScalar(rb.OutputType.segments, rb.ScalarType.radius))
    
    N = len(rs_segI)
    print("Number of segments is "+str(N))
    rs_segPot = np.ones(N) * -10000 # pressure head
    
    rs_segV = np.zeros((N,3)) # needed for doussan
    rs_segZ = np.zeros(N) # segment mid point z coordinate
    for i in range(0,N):
        s = rs_segI[i]
        n1 = v2v(rs_nodes[s.x])
        n2 = v2v(rs_nodes[s.y])
        rs_segV[i,:] = n2-n1
        mid = 0.5*(n1+n2)
        rs_segZ[i] = mid[2]  
    
    rs_flux = np.zeros(N)
    for i in range(0,N):
        rs_flux[i] = 0 # todo                  
     
    #
    # Apply Sink
    # 
     
              
    # Output
    print("time =", int(time), "\tdt =", dt, "\tIter. =", int(nrIterations),  "\tInf:", format(sumInfiltration, '.3f'))
    myPlot[0].clear()
    myPlot[0].set_xlim(0, 0.5)
    myPlot[0].set_xlabel("Water content [m$^3$ m$^{-3}$]",fontsize=16,labelpad=8)
    myPlot[0].set_ylabel("Depth [m]",fontsize=16,labelpad=8)
    myPlot[0].plot(inf.theta, -inf.z, 'k-')
    myPlot[0].plot(inf.theta, -inf.z, 'ko')
    myPlot[1].plot(time, flux, 'ko')
    plt.pause(0.0001)
            
    if (float(nrIterations/inf.maxNrIterations) < 0.1): 
        dt = min(dt*2, maxTimeStep) # go fasetr        
           
plt.ioff()
plt.show()
