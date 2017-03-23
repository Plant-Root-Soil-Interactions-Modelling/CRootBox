#PSP_vapor1D
from __future__ import division

import PSP_grid as grid
from PSP_ThomasAlgorithm import ThomasBoundaryCondition
from PSP_soil import *

waterDensity = 1000.        
area = 1                    
maxNrIterations = 100
tolerance = 1e-9

#vectors of numerical solution
n = 100                        
z = np.zeros(n+2, float)      
vol = np.zeros(n+2, float)     
a = np.zeros(n+2, float)      
b = np.zeros(n+2, float)       
c = np.zeros(n+2, float)      
d = np.zeros(n+2, float)      

dz = np.zeros(n+2, float)      
psi = np.zeros(n+2, float)    
dpsi = np.zeros(n+2, float)    
theta = np.zeros(n+2, float)   
oldTheta = np.zeros(n+2, float) 
vapor = np.zeros(n+2, float)  
oldvapor = np.zeros(n+2, float)  
C = np.zeros(n+2, float)      
k = np.zeros(n+2, float)       

u = np.zeros(n+2, float)       
du = np.zeros(n+2, float)      
f = np.zeros(n+2, float)       
H = np.zeros(n+2, float)       
H0 = np.zeros(n+2, float)     

def initializeWater(funcType, soil, theta_0):
    global z
    
    # vector depth [m]
    z = grid.linear(n, soil.lowerDepth)
    vol[0] = 0
    for i in range(n+1): 
        dz[i] = z[i+1]-z[i]
        if (i > 0): vol[i] = area * (z[i+1] - z[i-1]) / 2.0
    
    #initial conditions
    psi_0 = waterPotential(funcType, soil, theta_0)
    k_0 = (hydraulicConductivityFromTheta(funcType, soil, theta_0) 
            +  vaporConductivityFromPsiTheta(soil, psi_0, theta_0))

    psi[0] = 0
    for i in range(1, n+2):
        oldTheta[i] = theta_0
        theta[i] = theta_0
        oldvapor[i] = vaporFromPsi(soil, psi[i], theta[i])
        vapor[i] = vaporFromPsi(soil, psi[i], theta[i])
        psi[i] = psi_0
        H[i] = psi[i] - z[i]*g
        k[i] = k_0
            
def NewtonRapsonMP(funcType, soil, dt, isFreeDrainage):
    airEntry = airEntryPotential(funcType, soil)
    
    if (isFreeDrainage):
        psi[n+1] = psi[n]
        theta[n+1] = theta[n] 
        k[n+1] = k[n]
    
    nrIterations = 0
    massBalance = 1
    while ((massBalance > tolerance) and (nrIterations < maxNrIterations)):
        massBalance = 0
        for i in range(1, n+1):
            k[i] = (hydraulicConductivityFromTheta(funcType, soil, theta[i]) 
                    +  vaporConductivityFromPsiTheta(soil, psi[i], theta[i]))
            u[i] = g * k[i]
            du[i] = -u[i] * soil.Campbell_n / psi[i]
            capacity = dTheta_dPsi(funcType, soil, psi[i])
            capacity_vapor = dvapor_dPsi(funcType, soil, psi[i], theta[i])
            C[i] = (vol[i] * (waterDensity*capacity + capacity_vapor)) / dt
        
        for i in range (1, n+1):
            f[i] = ((psi[i+1] * k[i+1] - psi[i] * k[i]) 
                    / (dz[i] * (1 - soil.Campbell_n))) - u[i]
            if (i == 1):
                a[i] = 0.
                c[i] = -k[i+1] / dz[i]
                b[i] =  k[i] / dz[i] + C[i] + du[i]
                d[i] = evaporation_flux(psi[i]) - f[i] + vol[i]*(waterDensity 
                           * (theta[i] - oldTheta[i]) + (vapor[i]-oldvapor[i])) /dt
            else:
                a[i] = -k[i-1] / dz[i-1] - du[i-1]
                c[i] = -k[i+1] / dz[i]
                b[i] = k[i] / dz[i-1] + k[i] / dz[i] + C[i] + du[i]
                d[i] = f[i-1] - f[i] + (waterDensity * vol[i] 
                                        * (theta[i] - oldTheta[i])) /dt
            massBalance += abs(d[i])
    
        ThomasBoundaryCondition(a, b, c, d, dpsi, 1, n)
        
        for i in range(1, n+1):
            psi[i] -= dpsi[i] 
            psi[i] = min(airEntry, psi[i])    
            theta[i] = thetaFromPsi(funcType, soil, psi[i])
            vapor[i] = vaporFromPsi(soil, psi[i], theta[i])
        nrIterations += 1
        
        if (isFreeDrainage):
            psi[n+1] = psi[n]
            theta[n+1] = theta[n] 
            k[n+1] = k[n]
     
    if (massBalance < tolerance):
        flux = evaporation_flux(psi[1])
        return True, nrIterations, flux
    else:
        return False, nrIterations, 0