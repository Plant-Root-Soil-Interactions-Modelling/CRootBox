#PSP_vapor1D
from __future__ import print_function, division
from PSP_public import *
import PSP_grid as grid
import PSP_plant as plant
from PSP_ThomasAlgorithm import ThomasBoundaryCondition
from PSP_soil import *


z = np.zeros(n+2, float)      
vol = np.zeros(n+2, float)     
a = np.zeros(n+2, float)       
b = np.zeros(n+2, float)       
c = np.zeros(n+2, float)       
d = np.zeros(n+2, float)       

dz = np.zeros(n+2, float)           
psi = np.zeros(n+2, float)         
oldpsi = np.zeros(n+2, float)      
dpsi = np.zeros(n+2, float)        
theta = np.zeros(n+2, float)      
oldTheta = np.zeros(n+2, float)    
vapor = np.zeros(n+2, float)        
oldvapor = np.zeros(n+2, float)    
C = np.zeros(n+2, float)           
k = np.zeros(n+2, float)           
k_mean = np.zeros(n+2, float)      

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
           + vaporConductivityFromPsiTheta(soil, psi_0, theta_0))

    psi[0] = 0
    for i in range(1, n+2):
        oldTheta[i] = theta_0
        theta[i] = theta_0
        oldvapor[i] = vaporFromPsi(funcType, soil, psi[i], theta[i])
        vapor[i] = vaporFromPsi(funcType, soil, psi[i], theta[i])
        oldpsi[i] = psi_0
        psi[i] = psi_0
        H[i] = psi[i] - z[i]*g
        k[i] = k_0
    return
            
     
def NewtonRapsonMP(funcType, soil, dt, isFreeDrainage):
    if (isFreeDrainage):
        psi[n+1] = psi[n]
        theta[n+1] = theta[n] 
        k[n+1] = k[n]
    nrIterations = 0
    massBalance = 1.
    while ((massBalance > tolerance) and (nrIterations < maxNrIterations)):
        massBalance = 0
        for i in range(1, n+1):
            k[i] = (hydraulicConductivityFromTheta(funcType, soil, theta[i]) 
                +  vaporConductivityFromPsiTheta(soil, psi[i], theta[i]))
            k_mean[i] = kMean(LOGARITHMIC, k[i], k[i+1])
            u[i] = g * hydraulicConductivityFromTheta(funcType, soil, theta[i])
            du[i] = -u[i] * soil.Campbell_n / psi[i]
            capacity = dTheta_dPsi(funcType, soil, psi[i])
            capacity_vapor = dvapor_dPsi(funcType, soil, psi[i], theta[i])
            C[i] = (vol[i] * (waterDensity*capacity + capacity_vapor)) / dt
        for i in range (1, n+1):
            f[i] = (k_mean[i]*(psi[i+1]-psi[i])/dz[i]) - u[i]
            if (i == 1):
                a[i] = 0.
                c[i] = -k_mean[i] / dz[i]
                b[i] =  k_mean[i] / dz[i] + C[i] + du[i]
                d[i] = ((evaporation_flux(psi[i]) - f[i] + 
                vol[i]*(waterDensity * (theta[i] - oldTheta[i]) 
                + (vapor[i]-oldvapor[i])) /dt + plant.t[i]))
            else:
                a[i] = -k_mean[i-1] / dz[i-1] - du[i-1]
                c[i] = -k_mean[i] / dz[i]
                b[i] = k_mean[i-1] / dz[i-1] + k_mean[i] / dz[i] + C[i] + du[i]
                d[i] = ((f[i-1] - f[i] + 
                        vol[i]*(waterDensity * (theta[i] - oldTheta[i]) 
                           + (vapor[i]-oldvapor[i])) /dt + plant.t[i]))
                massBalance += abs(d[i])
        ThomasBoundaryCondition(a, b, c, d, dpsi, 1, n)
        for i in range(1, n+1):
            psi[i] -= dpsi[i]
            theta[i] = thetaFromPsi(funcType, soil, psi[i])
            vapor[i] = vaporFromPsi(funcType, soil, psi[i], theta[i])
        nrIterations += 1
        if (isFreeDrainage):
            psi[n+1] = psi[n]
            theta[n+1] = theta[n] 
            k[n+1] = k[n]
    if (massBalance < tolerance):
        flux = evaporation_flux(psi[1]) + sum(plant.t)
        return True, nrIterations, flux
    else:
        return False, nrIterations, 0
    
    
def cellCentFiniteVolWater(funcType, soil, dt, isFreeDrainage):
    if (isFreeDrainage):
        psi[n+1] = psi[n]
        theta[n+1] = theta[n] 
        k[n+1] = k[n]
    for i in range(1, n+1):
        theta[i] = oldTheta[i]
    massBalance = 1.
    nrIterations = 0
    while ((massBalance > tolerance) and (nrIterations < maxNrIterations)):
        for i in range(1, n+1):
            k[i] = (hydraulicConductivityFromTheta(funcType, soil, theta[i]) 
                   + vaporConductivityFromPsiTheta(soil, psi[i], theta[i]))
            capacity = dTheta_dPsi(funcType, soil, psi[i])
            capacity_vapor = dvapor_dPsi(funcType, soil, psi[i], theta[i])
            C[i] = (vol[i] * (waterDensity*capacity + capacity_vapor)) / dt
        f[0] = 0
        for i in range(1, n+1):
            f[i] = area * kMean(LOGARITHMIC, k[i], k[i+1]) / dz[i]
        for i in range(1, n+1):    
            a[i] = -f[i-1]
            if (i == 1):
                b[i] = C[i] + f[i]
                c[i] = -f[i]
                d[i] = (-evaporation_flux(psi[1]) + 
                        C[i] * oldpsi[i] - area * f[i]*dz[i] * g - plant.t[i])
            elif (i < n):
                b[i] = C[i] + f[i-1] + f[i]
                c[i] = -f[i]
                d[i] = (C[i] * oldpsi[i] 
                        - area * (f[i]*dz[i]-f[i-1]*dz[i-1]) * g - plant.t[i])
            elif (i == n):
                c[n] = 0
                if (isFreeDrainage):
                    b[n] = C[n] + f[n-1] 
                    d[n] = (C[n] * oldpsi[n] - 
                    area * (k[n] - f[i-1]*dz[i-1]) * g - plant.t[n])
                else:
                    b[n] = C[n] + f[n-1] + f[n]
                    d[n] = (C[n] * oldpsi[n] + f[n]*psi[n+1] 
                            - area * (f[i]*dz[i]-f[i-1]*dz[i-1]) * g - plant.t[i])
        ThomasBoundaryCondition(a, b, c, d, psi, 1, n)
        for i in range(1, n+1):
            theta[i] = thetaFromPsi(funcType, soil, psi[i])
            vapor[i] = vaporFromPsi(funcType, soil, psi[i], theta[i])
        if (isFreeDrainage):
            psi[n+1] = psi[n]
            theta[n+1] = theta[n] 
            k[n+1] = k[n]
        newSum = 0
        for i in range(1, n+1):
            newSum += C[i]*(psi[i]-oldpsi[i])
        if (isFreeDrainage):    
            massBalance = abs(newSum + evaporation_flux(psi[1]) + area*k[n]*g)
        else:
            massBalance = (abs(newSum + evaporation_flux(psi[1]) 
                        + f[n]*(psi[n]-psi[n+1]) + area*f[i]*dz[i]*g))
        nrIterations += 1
        
    if (massBalance < tolerance) or (dt == minTimeStep):
        flux = evaporation_flux(psi[1])+ sum(plant.t)
        return True, nrIterations, flux
    else:
        return False, nrIterations, 0
    
