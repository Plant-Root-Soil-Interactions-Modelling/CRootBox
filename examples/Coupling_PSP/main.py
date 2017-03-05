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
    
    
#
# look up function for soil matric potential
#    
def soil_p(x,y,z,psi,depth):
    i = round(-z/depth*inf.n)
    # print(i, z, inf.n, depth)
    return psi[i+1]
     
     
        
#
# Initialize soil domain (from Soil Physics with Python)
#    
isSuccess, soil = inf.readSoil("clay.txt")
if not isSuccess: 
    print("warning: wrong soil file.")
    quit()        
    
funcType = inf.CAMPBELL # inf.CAMPBELL, inf.RESTRICTED_VG, inf.IPPISCH_VG
        
Se = 0.5  # Initial degree of saturation ]0-1]

inf.initializeWater(funcType, soil, Se, inf.CELL_CENT_FIN_VOL) # for simplicity we use a linear grid
        
# ubPotential = inf.airEntryPotential(funcType, soil[0]) # [J kg^-1] upper boundary condition INFILTRATION
ubPotential = inf.waterPotential(funcType, soil[0],inf.thetaFromSe(funcType, soil[0],0.3))    
print("Upper boundary potential ", ubPotential)
isFreeDrainage = True

sumInfiltration = 0
totalIterationNr = 0

#
# Initialize root domain
#
# rsname = "Anagallis_femina_Leitner_2010" 
rsname = "Anagallis_femina_Leitner_2010"
rs = rb.RootSystem()
rs.openFile(rsname,parameterPath())
rs.initialize() # hydrotropism is not set right now, link to soil is missing

rs_Kr = np.array([ 1.16e-6, 1.74e-5, 1.74e-5, 1.74e-5, 11.74e-5, 1.74e-5, 1.74e-5 ]) # s/m; root hydraulic radial conductivity per root type 
rs_Kz = np.array([ 2.3e-8, 1.16e-11, 1.16e-11, 1.16e-11, 1.16e-11, 1.16e-11, 1.16e-11 ]) # m²*s; root hydraulic axial conductivity per root type  

dirichlet = False

#
# Simulation parameters
#                   
simTime = 30*24*3600   
maxTimeStep = 360           
dt = 36              
time = 0               

out_delay = 24*3600 # figure output delay
out_next = out_delay             
out_img = np.zeros((inf.n+2,round(simTime/out_delay)))
out_c = 0

cflux = 0 # cumulative TODO
  
#   
# Prepare output
#
plt.ion()
fig = plt.figure() 

ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2, projection='3d')  
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)

#
# Simulation loop
#
while (time < simTime):
    
    dt = min(dt, out_next - time) 
    print()
    print("Time step "+ str(dt))

    #
    # Richards Code (from Soil Physics with Python)
    #     
    t1 = timer.time()
    print("1. Richards equation")   
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
    t2 = timer.time();
    print("2 . Root growth")
    rs.simulate(dt/3600/24) # seconds to days
        
    #
    # Root System fluxes
    #    
    t3 = timer.time()
    print("3 . Xylem fluxes")
    seg = seg2a(rs.getSegments())
    nodes = vv2a(rs.getNodes())/100. # convert to meter
    rs_ana = rb.SegmentAnalyser(rs) 
    type = v2a(rs_ana.getScalar(rb.ScalarType.type))
    radius = v2a(rs_ana.getScalar(rb.ScalarType.radius))/100. # convert to meter 
    time_ = v2a(rs_ana.getScalar(rb.ScalarType.time))*3600*24 # convert to seconds 
    
    kr = np.array(list(map(lambda t: rs_Kr[int(t)-1], type))) # convert from 'per type' to 'per segment'
    kz = np.array(list(map(lambda t: rs_Kz[int(t)-1], type)))
    # age dependency
    #kr = kr * (10*time_ + 1) 
    #kz = kz / (10*time_ + 1)
    
    rho = 1e3 # kg / m^3      
    g = 9.8  # m / s^2    
     
    soil_p2 = lambda x,y,z : soil_p(x,y,z, inf.psi,soil[-1].lowerDepth) # J/kg
    
    Q, b = xylem_flux.linear_system(seg, nodes, radius, kr, kz, rho, g, soil_p2)  
    
    Tpot = np.array([-5.79e-4]) # m^3 s^-1 

    
    if dirichlet == False:
        Q, b = xylem_flux.bc_neumann(Q, b, np.array([0]), Tpot, seg, nodes)
        x = LA.spsolve(Q, b) # direct
        Teff = xylem_flux.axial_flux0(x, seg, nodes, kz, rho, g) # verify that Teff == Tpot        
        print("using neumann ( Effective Transpiration = " +str(Teff)+" m^3 s^-1)" )        
        
    if dirichlet or (x[0]<-1500):        
        dirichlet = True
        Q, b = xylem_flux.bc_dirichlet(Q, b, np.array([0]), np.array([-1500])) # J/kg
        x = LA.spsolve(Q, b) # direct
        Teff = xylem_flux.axial_flux0(x, seg, nodes, kz, rho, g)
        dirichlet = Teff>Tpot
        print("using dirichlet ( Effective Transpiration = "+str(Teff)+" m^3 s^-1)")
        
    radial_flux = xylem_flux.radial_flux(x, seg, nodes, radius, kr, soil_p2)    
    
    print("Minimal xylem potential", min(x), " Maximal xylem potential ", max(x))
    
    #
    # Sink term
    #    
    t4 = timer.time()
    print("4 . Apply sink")
    sink = np.zeros(inf.n)
    for i in range(0,len(seg)):
        z1 = nodes[seg[i,0],2]
        z2 = nodes[seg[i,1],2]
        z = 0.5*(z1+z2)
        ind = round(-z/soil[-1].lowerDepth*inf.n)
        sink[ind] += radial_flux[i] * dt       
    
    for i in range(0,inf.n):
        inf.theta[i+1] += sink[i] 
        inf.theta[i+1] = min(inf.theta[i+1],1)
        inf.theta[i+1] = max(inf.theta[i+1],0)          
      
    # update psi
    for i in range(0,inf.n+2):
        inf.psi[i] = inf.waterPotential(funcType, soil[inf.hor[i]], inf.theta[i])
     
    cflux += -float(np.sum(sink))   
     
    #
    # Output
    #
    t = timer.time()
    #print("flux "+ str(np.sum(sink)))
    t0 = t-t1    
    s0 = str(format(t0, '.3f'))
    s1 = str(format((t2-t1)/t0, '.3f'))
    s2 = str(format((t3-t2)/t0, '.3f'))
    s3 = str(format((t4-t3)/t0, '.3f'))
    s4 = str(format((t -t4)/t0, '.3f'))
    print("simtime =", time/3600./24. , "days, spent = "+ s0 +"s ("+s1+", "+s2+", "+s3+", "+s4+") ", rs.getNumberOfNodes(), "nodes, ", " dt =", dt, " iterations =", int(nrIterations),  " summed infiltration:", format(sumInfiltration, '.3f'))
    
    if time>=out_next:
        out_next += out_delay
        
        # Subplot 1 - WATER CONTENT
        out_img[:,out_c] = inf.theta[:]
        imin = np.min(np.min(out_img))
        imax = np.max(np.max(out_img))     
        out_c += 1
        #ax1.imshow((out_img-imin)/(imax-imin), interpolation="nearest", cmap="hot")        
        ax1.set_xlabel("Time");     
        ax1.set_ylabel("Depth");
        ax1.set_title("Water content")
        # ax1.yticks(np.linspace(0,-1,11)) # how?        
        ax1.imshow(out_img, interpolation="nearest", cmap="hot")     
        # fig.colorbar(ax1,[-1,0,1]) # np.linspace(imin,imax,10) 

#        Subplot 4 - THETA        
#         ax1.set_ylim(-1, 0)
#         ax1.set_xlabel("Water potential");
#         ax1.set_ylabel("Depth [m]")         
#         ax1.plot(inf.psi,-inf.z)   
        
        # Subplot2 - ROOT SYSTEM
        plotRSscatter(ax2, rs.getRootTips())        
        
        # Subplot3  -SINK
        ax3.set_ylim(-1, 0)
        ax3.set_xlabel("Sink");
        ax3.set_ylabel("Depth [m]")         
        ax3.plot(sink,-inf.z[1:len(inf.z)-1])         
        
#         # Subplot 3 - CUMULATIVE ROOT WATER UPTAKE        
#         ax3.set_xlabel("Time [days]")
#         ax3.set_ylabel("Root water uptake [1]") 
#         ax3.plot([time/3600/24,out_next/3600/24], [cflux,cflux], 'k-')
    
        # Subplot 4 - THETA
        ax4.set_ylim(-1, 0)
        ax4.set_xlabel("Water content");
        ax4.set_ylabel("Depth [m]")         
        ax4.plot(inf.theta,-inf.z)   
    
        plt.pause(0.0001)
            
    if (float(nrIterations/inf.maxNrIterations) < 0.1): 
        dt = min(dt*2, maxTimeStep) # go fasetr        
           
plt.ioff()
plt.show()
