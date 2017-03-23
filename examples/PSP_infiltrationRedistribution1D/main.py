#PSP_infiltrationRedistribution1D
from __future__ import print_function, division

import matplotlib.pyplot as plt
import PSP_infiltration1D as inf
    
def main():  
    isSuccess, soil = inf.readSoil("siltLoam.txt")
    if not isSuccess: 
        print("warning: wrong soil file.")
        return
    
#     print (inf.CAMPBELL,' Campbell')
#     print (inf.RESTRICTED_VG,' van Genuchten with m = 1-1/n restriction')
#     print (inf.IPPISCH_VG,' Ippisch-van Genuchten')
#     funcType = int(input("Select water retention curve: "))
    funcType = 1    
    
#     if (funcType == inf.CAMPBELL) and (len(soil) == 1):
#         print()
#         print (inf.CELL_CENT_FIN_VOL,' Cell-Centered Finite Volume')
#         print (inf.NEWTON_RAPHSON_MP,' Matric Potential with Newton-Raphson')
#         print (inf.NEWTON_RAPHSON_MFP,' Matric Flux Potential with Newton-Raphson')
#         solver = int(input("Select solver: "))
#     else:
#         solver = inf.CELL_CENT_FIN_VOL
    solver = inf.CELL_CENT_FIN_VOL    
    
#     myStr = "Initial degree of saturation ]0-1]:" 
#     Se = inf.NODATA
#     print()
#     while ((Se <= 0.0) or (Se > 1.0)):
#         Se = float(input(myStr))
    Se = 0.4

#     print()
#     print ("1: Free drainage")
#     print ("2: Constant water potential")
#     boundary = int(input("Select lower boundary condition:"))
#     if (boundary == 1):
#         isFreeDrainage = True
#     else:
#         isFreeDrainage = False
    isFreeDrainage = True

    inf.initializeWater(funcType, soil, Se, solver)
            
    ubPotential = inf.airEntryPotential(funcType, soil[0])  # [J kg^-1] upper boundary condition    
                
     
    # simulationLenght = int(input("\nNr of simulation hours:")) # hours of simulation        
    simulationLenght = 24*30
                    
    endTime = simulationLenght * 3600   
    maxTimeStep = 3600                  
    dt = maxTimeStep / 10               
    time = 0                            
    sumInfiltration = 0
    totalIterationNr = 0
         
    plt.ion()
    f, myPlot = plt.subplots(2, figsize=(10, 8), dpi=80)
    myPlot[1].set_xlim(0, simulationLenght * 3600)
    myPlot[1].set_ylim(0, 0.1)
    myPlot[1].set_xlabel("Time [s]",fontsize=16,labelpad=8)
    myPlot[1].set_ylabel("Infiltration Rate [kg m$^{-2}$ s$^{-1}$]",fontsize=16,labelpad=8)
    plt.tick_params(axis='both', which='major', labelsize=12,pad=6)

    while (time < endTime):
        dt = min(dt, endTime - time)
        if (solver == inf.CELL_CENT_FIN_VOL):
            success, nrIterations, flux = inf.cellCentFiniteVolWater(funcType, 
                        soil, dt, ubPotential, isFreeDrainage,inf.LOGARITHMIC)
        elif (solver == inf.NEWTON_RAPHSON_MP):
            success, nrIterations, flux = inf.NewtonRapsonMP(funcType,
                        soil, dt, ubPotential, isFreeDrainage)
        elif (solver == inf.NEWTON_RAPHSON_MFP):
            success, nrIterations, flux = inf.NewtonRapsonMFP(funcType,
                        soil, dt, ubPotential, isFreeDrainage)
        totalIterationNr += nrIterations
        
        if success:
            for i in range(inf.n+2):
                inf.oldTheta[i] = inf.theta[i]
            sumInfiltration += flux * dt 
            time += dt
                        
            print("time =", int(time), "\tdt =", dt, 
                  "\tIter. =", int(nrIterations), 
                  "\tInf:", format(sumInfiltration, '.3f'))
            myPlot[0].clear()
            myPlot[0].set_xlim(0, 0.5)
            myPlot[0].set_xlabel("Water content [m$^3$ m$^{-3}$]",fontsize=16,labelpad=8)
            myPlot[0].set_ylabel("Depth [m]",fontsize=16,labelpad=8)
            myPlot[0].plot(inf.theta, -inf.z, 'k-')
            myPlot[0].plot(inf.theta, -inf.z, 'ko')
            myPlot[1].plot(time, flux, 'ko')
            plt.pause(0.0001)
            
            if (float(nrIterations/inf.maxNrIterations) < 0.1): 
                    dt = min(dt*2, maxTimeStep)
        
        else:
            print ("dt =", dt, "No convergence")
            dt = max(dt / 2, 1)
            for i in range(inf.n+2):
                inf.theta[i] = inf.oldTheta[i]
                if (solver == inf.NEWTON_RAPHSON_MFP):
                    inf.psi[i] = inf.MFPFromTheta(soil[inf.hor[i]], inf.theta[i])
                else:
                    inf.psi[i] = inf.waterPotential(funcType, soil[inf.hor[i]], inf.theta[i])
           
    print("nr of iterations per hour:", totalIterationNr / simulationLenght)
    plt.ioff()
    plt.show()
main()
