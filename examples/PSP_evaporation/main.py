#PSP_evaporation
from __future__ import print_function, division
import matplotlib.pyplot as plt
import PSP_vapor1D as vap
    
def main():  
    isSuccess, soil = vap.readSoil("siltLoam.txt")
    if not isSuccess: 
        print("warning: wrong soil file.")
        return
    
    funcType = vap.CAMPBELL
    
    #Initial conditions 
    thetaIni = 0.2         
    isFreeDrainage = True
    
    simulationLenght = 80  
    vap.initializeWater(funcType, soil, thetaIni)
    endTime = simulationLenght * 3600   
    maxTimeStep = 3600                  
    dt = 600                           
    time = 0                           
    sumEvaporation = 0
    
    plt.ion()
    f, myPlot = plt.subplots(2, figsize=(8, 8), dpi=80)
    myPlot[0].set_xlim(0, soil.thetaS)
            
    myPlot[1].set_xlim(0, simulationLenght)
    myPlot[1].set_ylim(0, 0.3)
    myPlot[1].set_ylabel("Evaporation Rate [mm h$^{-1}$]",fontsize=16,labelpad=8)
    myPlot[1].set_xlabel("Time [h]",fontsize=16,labelpad=8)
    
    while (time < endTime):
        dt = min(dt, endTime - time)
        success, iterations, evaporationFlux = vap.NewtonRapsonMP(funcType, soil, dt, isFreeDrainage)
        
        if success:
            for i in range(vap.n+1):
                vap.oldTheta[i] = vap.theta[i]
                vap.oldvapor[i] = vap.vapor[i]
            sumEvaporation += evaporationFlux * dt 
            time += dt
                        
            print("time =", int(time), "\tdt =", int(dt), "\tIter. =", iterations, 
                  "\tsumEvap:", format(sumEvaporation, '.3f'))
            
            myPlot[0].clear()
            myPlot[0].set_xlim(0, soil.thetaS)
            myPlot[0].set_xlabel("Water content [m$^3$ m$^{-3}$]",fontsize=16,labelpad=8)
            myPlot[0].set_ylabel("Depth [m]",fontsize=16,labelpad=8)
            myPlot[0].plot(vap.theta, -vap.z, 'yo')
            myPlot[1].plot(time/3600., evaporationFlux*3600., 'ro')
            plt.pause(0.0001)
            
            if (float(iterations / vap.maxNrIterations) < 0.1): 
                    dt = min(dt*2.0, maxTimeStep)  
        else:
            print("time =", int(time), "\tdt =", int(dt), 
                  "\tIter. =", iterations, "No convergence")
            
            for i in range(vap.n+1):
                vap.theta[i] = vap.oldTheta[i]
                vap.vapor[i] = vap.oldvapor[i]
                vap.psi[i] = vap.waterPotential(funcType, soil, vap.theta[i])
            dt = max(dt / 2.0, 1)
    
    plt.ioff()        
    plt.show()
main()
