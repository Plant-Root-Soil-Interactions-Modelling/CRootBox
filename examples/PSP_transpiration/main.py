#PSP_transpiration
from __future__ import print_function, division
from PSP_public import *
import matplotlib.pyplot as plt
import PSP_vapor1D as vap

def main():  
    isSuccess, soil = vap.readSoil("soil.txt")
    if not isSuccess: 
        print("warning: wrong soil file.")
        return
    
    funcType = vap.CAMPBELL
    
#     print (vap.CELL_CENT_FIN_VOL,' Cell-Centered Finite Volume')
#     print (vap.NEWTON_RAPHSON_MP,' Matric Potential with Newton-Raphson')
#     solver = int(input("Select solver: "))
    solver =vap.CELL_CENT_FIN_VOL

    myStr = "]0, " + format(soil.thetaS, '.2f') 
    myStr += "] initial water content (m^3 m^-3):" 
    thetaIni = vap.NODATA
    thetaIni = 0.3
    print("thetaIni", thetaIni, soil.VG_thetaR, soil.thetaS)
#     print()
#     while ((thetaIni <= soil.VG_thetaR) or (thetaIni > soil.thetaS)):
#         thetaIni = float(input(myStr))

    vap.initializeWater(funcType, soil, thetaIni)
    vap.plant.InitTransp(soil, vap.z)  
  
#     print()
#     print ("1: Free drainage")
#     print ("2: Constant water potential")
#     boundary = int(input("Select lower boundary condition:"))
#     if (boundary == 1):
#         isFreeDrainage = True
#     else:
#         isFreeDrainage = False
    isFreeDrainage = True
   
#     simulationLenght = int(input("\nNr of simulation hours:"))    
    simulationLenght = 14*24
    endTime = simulationLenght * 3600 # in seconds
    dt = maxTimeStep / 10               
    time = 0                            
    sumETr = 0
    totalIterationNr = 0
	
    plt.ion()
    f, myPlot = plt.subplots(4, figsize=(8, 9), dpi=80)
    f.subplots_adjust(hspace=.35)
    myPlot[0].tick_params(axis='both', which='major', labelsize=12,pad=4)
    myPlot[1].tick_params(axis='both', which='major', labelsize=12,pad=4)
    myPlot[2].tick_params(axis='both', which='major', labelsize=12,pad=4)
    myPlot[3].tick_params(axis='both', which='major', labelsize=12,pad=4)
    
    myPlot[1].set_xlim(0, simulationLenght)
    myPlot[1].set_ylim(0, 0.8)
    myPlot[1].set_ylabel("Evapotransp. Rate [mm h$^{-1}$]")
    myPlot[1].set_xlabel("Time [h]")
    myPlot[3].set_xlim(0, simulationLenght)
    myPlot[3].set_xlabel("Time [-]")
    myPlot[3].set_ylabel("Water potential [J kg$^{-1}$]")
    
    while (time < endTime):
    
        dt = min(dt, endTime - time)
        TimeOfDay = float(time%(3600*24))/3600.
        tp = vap.fi * vap.ETp * 2.3 * vap.np.power(0.05 + vap.np.sin(0.0175 * 7.5 * TimeOfDay),4.)           
        vap.plant.PlantWaterUptake(tp, funcType, soil, vap.theta, vap.psi, vap.z)
        
        if (solver == vap.CELL_CENT_FIN_VOL):
            success, nrIterations, flux = vap.cellCentFiniteVolWater(funcType, soil, dt, isFreeDrainage)
        elif (solver == vap.NEWTON_RAPHSON_MP):
            success, nrIterations, flux = vap.NewtonRapsonMP(funcType, soil, dt, isFreeDrainage)
        totalIterationNr += nrIterations
        
        transp, leafPot, soilPot=vap.plant.PlantWaterUptake(tp, funcType, soil, vap.theta, vap.psi, vap.z) 
        print ("leafPot",leafPot)
        print ("transp",transp)        
        
        if (success):
            for i in range(vap.n+1):
                vap.oldTheta[i] = vap.theta[i]
                vap.oldvapor[i] = vap.vapor[i]
                vap.oldpsi[i] = vap.psi[i]
                
            sumETr += flux * dt 
            print("time =", int(time), "\tdt =", int(dt), 
                  "\tIter. =", int(nrIterations), 
                  "\tsum ETr:", format(sumETr, '.3f'))
            time += dt
          
            myPlot[0].clear()
            myPlot[0].set_xlabel("Water content [m$^3$ m$^{-3}$]")
            myPlot[0].set_ylabel("Depth [m]")
            myPlot[0].set_xlim(-0, soil.thetaS)
            myPlot[0].plot(vap.theta[1:vap.n+1], -vap.z[1:vap.n+1], 'ko')
            myPlot[1].plot(time/3600., flux*3600., 'ko')
            myPlot[2].clear()
            myPlot[2].set_xlabel("Root Water Uptake [-]")
            myPlot[2].set_ylabel("Depth [m]")
            myPlot[2].set_xlim(-0.5e-5,1.5e-5)
            myPlot[2].plot(vap.plant.t[1:vap.n+1], -vap.z[1:vap.n+1], 'ko')
            myPlot[3].plot(time/3600., leafPot, 'ko')
            myPlot[3].plot(time/3600., soilPot, 'k^')
            plt.pause(0.0001)
            
            if (float(nrIterations/vap.maxNrIterations) < 0.1): 
                    dt = min(dt*2, maxTimeStep)
        else:
            print ("time =", int(time), "\tdt =", int(dt), "No convergence")
            dt = max(dt / 2, minTimeStep)
            for i in range(vap.n+1):
                vap.theta[i] = vap.oldTheta[i]
                vap.vapor[i] = vap.oldvapor[i]
                vap.psi[i] = vap.waterPotential(funcType, soil, vap.theta[i])
            
    plt.ioff()
    plt.show()
main()

