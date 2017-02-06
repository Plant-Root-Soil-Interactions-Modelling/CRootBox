import math
import time as timer
import os # add the search path for py_rootbox.so (probably there is a nicer way to do it?)
import sys
cwd = os.getcwd()
i = cwd.index("CRootBox")
sys.path.append(cwd[0:i+8])

import numpy as np
import scipy.sparse.linalg as LA
from scipy import sparse

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import PSP_infiltration1D as inf

import py_rootbox as rb    
from rb_tools import *
import xylem_flux 
    
def soil_p(x,y,z,psi,depth):
    i = round(-z/depth*inf.n)
    return psi[i]
        
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
rsname = "anagallis_Leitner_et_al(2010)" 
rs = rb.RootSystem()
rs.openFile(rsname,parameterPath())
rs.initialize() # hydrotropism is not set right now, link to soil is missing

rs_Kr = ( 1e-10, 1e-10, 1e-10, 1e-10, 1e-10 ) # root hydraulic conductivity per root type 


#
# Simulation
#                   
simTime = 30*24*3600   
maxTimeStep = 24*3600                  
dt = 36              
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
ax1.set_xlabel("Time");     
ax1.set_ylabel("Depth");
ax1.set_title("Water content")
# ax1.yticks(np.linspace(0,-1,11)) # how?

ax2 = fig.add_subplot(2, 2, 2, projection='3d')  

ax3 = fig.add_subplot(2, 2, 3)
ax3.set_xlabel("Time [days]");
ax3.set_ylabel("Root water uptake [1]") 

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
    
    success, nrIterations, flux = inf.cellCentFiniteVolWater(funcType, soil, dt, ubPotential, isFreeDrainage, inf.LOGARITHMIC)
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
    t = timer.time()
    seg = seg2a(rs.getSegments())
    nodes = vv2a(rs.getNodes()) /100 # convert to meter
    rs_ana = rb.SegmentAnalyser(rs) 
    type = v2a(rs_ana.getScalar(rb.ScalarType.type))
    radius = v2a(rs_ana.getScalar(rb.ScalarType.radius))/100 # convert to meter 
    time_ = v2a(rs_ana.getScalar(rb.ScalarType.time))*3600*24 # convert to seconds 
    
    rs_Kr = np.array([ 1.16e-6, 1.74e-5, 1.74e-5, 1.74e-5, 11.74e-5, 1.74e-5, 1.74e-5 ]) # s/m; root hydraulic radial conductivity per root type 
    rs_Kz = np.array([ 2.3e-8, 1.16e-11, 1.16e-11, 1.16e-11, 1.16e-11, 1.16e-11, 1.16e-11 ]) # m²*s; root hydraulic axial conductivity per root type  
    kr = np.array(list(map(lambda t: rs_Kr[int(t)-1], type))) # convert from 'per type' to 'per segment'
    kz = np.array(list(map(lambda t: rs_Kz[int(t)-1], type)))
    #kr = kr * (10*time_ + 1)
    #kz = kz / (10*time_ + 1)
    
    rho = 1e3 # kg / m^3      
    g = 9.8  # m / s^2    
     
    soil_p2 = lambda x,y,z : soil_p(x,y,z, inf.psi,soil[-1].lowerDepth) # J/kg
    
    Q, b = xylem_flux.linear_system(seg, nodes, radius, kr, kz, rho, g, soil_p2)  
    Tpot = np.array([-5.79e-4]) # m^3 s^-1 
    Q, b = xylem_flux.bc_neumann(Q, b, np.array([0]), Tpot, seg, nodes)
    x = LA.spsolve(Q, b) # direct
    if x[0]<-1500:
        print("using dirichlet")
        Q, b = xylem_flux.bc_dirichlet(Q, b, np.array([0]), np.array([-1500])) # J/kg
        x = LA.spsolve(Q, b) # direct
        
    radial_flux = xylem_flux.radial_flux(x, seg, nodes, radius, kr, soil_p2)
    
    print("xylem flux time: " + str(timer.time()-t) +" sec" )

    #
    # Apply Sink
    #    
    sink = np.zeros(inf.n)
    for i in range(0,len(seg)):
        n1 = nodes[seg[i,0],:]
        n2 = nodes[seg[i,1],:]
        z = float(0.5*(n1[2]+n2[2]))
        ind = round(-z/soil[-1].lowerDepth*inf.n)
        sink[ind]+=radial_flux[i] *dt       
    
    for i in range(0,inf.n):
        inf.theta[i+1] -= sink[i] # TODO  
        inf.theta[i+1] = min(inf.theta[i+1],1)
        inf.theta[i+1] = max(inf.theta[i+1],0)          
      
    #print(inf.theta)
    for i in range(0,inf.n+2):
        inf.psi[i] = inf.waterPotential(funcType, soil[inf.hor[i]], inf.theta[i])
     
    print("flux "+ str(np.sum(sink)))
    cflux += -float(np.sum(sink))   
     
    #
    # Output
    #
    print("time =", int(time), "\tdt =", dt, "\tIter. =", int(nrIterations),  "\tInf:", format(sumInfiltration, '.3f'))
    
    if time>=out_next:
        out_next += out_delay
        
        # Subplot 1
        out_img[:,out_c] = inf.theta[:]
        imin = np.min(np.min(out_img))
        imax = np.max(np.max(out_img))     
        out_c += 1
        #ax1.imshow((out_img-imin)/(imax-imin), interpolation="nearest", cmap="hot")        
        ax1.imshow(out_img, interpolation="nearest", cmap="hot")     
        # fig.colorbar(ax1,[-1,0,1]) # np.linspace(imin,imax,10) 
        
        # Subplot2
        plotRSscatter(ax2, rs.getRootTips())        
        
        # Subplot 3                    
        ax3.plot([time/3600/24,out_next/3600/24], [cflux,cflux], 'k-')
    
        # Subplot 4
        ax4.plot(inf.theta,-inf.z)   
    
        plt.pause(0.0001)
            
    if (float(nrIterations/inf.maxNrIterations) < 0.1): 
        dt = min(dt*2, maxTimeStep) # go fasetr        
           
plt.ioff()
plt.show()
