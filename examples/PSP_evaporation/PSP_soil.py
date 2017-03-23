#PSP_soil.py
from __future__ import division

from PSP_readDataFile import *
from math import sqrt, log, exp

g = 9.8065                 
waterDensity = 1000.        
area = 1                    
mw = 0.018                 
R = 8.31                    
T = 293                     
dv = 0.000024              
vp = 0.017                  
E_p = 0.1*0.21/3600.            
h_a = 0.5                   
 
NODATA = -9999.

CAMPBELL = 1
RESTRICTED_VG = 2
IPPISCH_VG = 3
VAN_GENUCHTEN = 4

CELL_CENT_FIN_VOL = 1
NEWTON_RAPHSON_MP = 2
NEWTON_RAPHSON_MFP = 3

LOGARITHMIC = 0
HARMONIC = 1
GEOMETRIC = 2

class Csoil:
    upperDepth = NODATA       
    lowerDepth = NODATA       
    Campbell_he = NODATA       
    Campbell_b = NODATA        
    CampbellMFP_he = NODATA    
    VG_alpha = NODATA          
    VG_n = NODATA              
    VG_m = NODATA
    VG_he = NODATA             
    VG_alpha_mod = NODATA      
    VG_n_mod = NODATA          
    VG_m_mod = NODATA
    VG_Sc = NODATA            
    VG_thetaR = NODATA         
    Mualem_L = NODATA          
    thetaS = NODATA            
    Ks = NODATA                
    
def readSoil(soilFileName):
    soil = Csoil()
    A, isFileOk = readDataFile(soilFileName, 1, ',', False)
    if ((not isFileOk) or (len(A[0]) < 12)):
        return False, soil
    
    soil.upperDepth = A[0,0]
    soil.lowerDepth = A[0,1]
    soil.Campbell_he = A[0,2]
    soil.Campbell_b = A[0,3]
    soil.Campbell_n = 2.0 + (3.0 / soil.Campbell_b)
    soil.VG_he = A[0,4]
    soil.VG_alpha = A[0,5]
    soil.VG_n = A[0,6]
    soil.VG_m =  1. - (1. / soil.VG_n)
    soil.VG_alpha_mod = A[0,7]
    soil.VG_n_mod = A[0,8]
    soil.VG_m_mod =  1. - (1. / soil.VG_n_mod)
    soil.VG_Sc = (1. + (soil.VG_alpha_mod * abs(soil.VG_he))**soil.VG_n_mod)**(-soil.VG_m_mod)
    soil.VG_thetaR = A[0,9]
    soil.thetaS = A[0,10]
    soil.Ks = A[0,11]
    soil.Mualem_L = 0.5
    soil.CampbellMFP_he = soil.Ks * soil.Campbell_he / (1.0 - soil.Campbell_n) 
    return True, soil

def airEntryPotential(funcType, soil): 
    if (funcType == CAMPBELL):
        return(soil.Campbell_he)
    elif (funcType == IPPISCH_VG):
        return(soil.VG_he)
    elif (funcType == RESTRICTED_VG):
        return(0)
    else:
        return(NODATA)
     
def waterPotential(funcType, soil, theta):
    psi = NODATA
    Se = SeFromTheta(funcType, soil, theta)
    if (funcType == RESTRICTED_VG):
        psi = -(1./soil.VG_alpha)*((1./Se)**(1./soil.VG_m) - 1.)**(1./soil.VG_n)
    elif (funcType == IPPISCH_VG):
        psi = -(1./soil.VG_alpha_mod)*((1./(Se*soil.VG_Sc))**(1./soil.VG_m_mod)-1.)**(1./soil.VG_n_mod)
    elif (funcType == CAMPBELL):
        psi = soil.Campbell_he * Se**(-soil.Campbell_b)
    return(psi)
    
def SeFromTheta(funcType, soil, theta):
    if (theta >= soil.thetaS): return(1.)
    if (funcType == CAMPBELL):
        Se = theta / soil.thetaS
    else:
        Se = (theta - soil.VG_thetaR) / (soil.thetaS - soil.VG_thetaR)
    return (Se)

def thetaFromSe(funcType, soil, Se):
    if (funcType == RESTRICTED_VG) or (funcType == IPPISCH_VG):
        theta = (Se * (soil.thetaS - soil.VG_thetaR) + soil.VG_thetaR)
    elif (funcType == CAMPBELL):
        return(Se * soil.thetaS) 
    return(theta)

def degreeOfSaturation(funcType, soil, psi):
    if (psi >= 0.): return(1.)
    Se = NODATA
    if (funcType == IPPISCH_VG):
        if (psi >= soil.VG_he): Se = 1.
        else: 
            Se = (1./soil.VG_Sc) * pow(1.+pow(soil.VG_alpha_mod 
                                            * abs(psi), soil.VG_n_mod), -soil.VG_m_mod)
    elif (funcType == RESTRICTED_VG):
        Se = 1 / pow(1 + pow(soil.VG_alpha * abs(psi), soil.VG_n), soil.VG_m)  
    elif (funcType == CAMPBELL):
        if psi >= soil.Campbell_he: Se = 1.
        else: Se = pow(psi / soil.Campbell_he, -1. / soil.Campbell_b)
    return(Se) 

def thetaFromPsi(funcType, soil, psi):
    Se = degreeOfSaturation(funcType, soil, psi)
    theta = thetaFromSe(funcType, soil, Se)
    return(theta)
          
def hydraulicConductivityFromTheta(funcType, soil, theta): 
    k = NODATA      
    if (funcType == RESTRICTED_VG):
        Se = SeFromTheta(funcType, soil, theta)
        k = soil.Ks * pow(Se, soil.Mualem_L) * (1. -pow(1. -pow(Se, 1./soil.VG_m), soil.VG_m))**2
    elif (funcType == IPPISCH_VG):
        Se = SeFromTheta(funcType, soil, theta)
        num   = 1. - pow(1. - pow(Se * soil.VG_Sc, 1./ soil.VG_m_mod), soil.VG_m_mod);
        denom = 1. - pow(1. - pow(soil.VG_Sc, 1./ soil.VG_m_mod), soil.VG_m_mod);
        k = soil.Ks * pow(Se, soil.Mualem_L) * pow((num / denom), 2.)
    elif (funcType == CAMPBELL):
        psi = waterPotential(funcType, soil, theta)
        k = soil.Ks * (soil.Campbell_he / psi)**soil.Campbell_n
    return(k)

#---------------------------------------------
# dTheta/dH = dSe/dH (Theta_s - Theta_r)
#---------------------------------------------
def dTheta_dPsi(funcType, soil, psi):
    airEntry = airEntryPotential(funcType, soil)
    if (psi > airEntry): return 0.0
     
    if (funcType == RESTRICTED_VG):
        dSe_dpsi = soil.VG_alpha * soil.VG_n * (soil.VG_m 
                * pow(1. + pow(soil.VG_alpha * abs(psi), soil.VG_n), 
                -(soil.VG_m + 1.)) * pow(soil.VG_alpha * abs(psi), soil.VG_n - 1.))      
        return dSe_dpsi * (soil.thetaS - soil.VG_thetaR)
    elif (funcType == IPPISCH_VG):
        dSe_dpsi = soil.VG_alpha_mod * soil.VG_n_mod * (soil.VG_m_mod 
                * pow(1. + pow(soil.VG_alpha_mod * abs(psi), soil.VG_n_mod), 
                -(soil.VG_m_mod + 1.)) * pow(soil.VG_alpha_mod * abs(psi), soil.VG_n_mod - 1.))      
        dSe_dpsi *= (1. / soil.VG_Sc)
        return dSe_dpsi * (soil.thetaS - soil.VG_thetaR)
    elif (funcType == CAMPBELL):
        theta = soil.thetaS * degreeOfSaturation(funcType, soil, psi) 
        return -theta / (soil.Campbell_b * psi)

def MFPFromTheta(soil, theta): 
    return (soil.CampbellMFP_he * (theta / soil.thetaS)**(soil.Campbell_b + 3.0)) 

def MFPFromPsi(soil, psi):
    return (soil.CampbellMFP_he * (psi / soil.Campbell_he)**(1.0 - soil.Campbell_n)) 

def thetaFromMFP(soil, MFP):
    if (MFP > soil.CampbellMFP_he):
        return(soil.thetaS) 
    else:
        return(soil.thetaS * (MFP / soil.CampbellMFP_he)**(1.0/(soil.Campbell_b + 3.0)))  

def hydraulicConductivityFromMFP(soil, MFP):
    b3 = (2.0 * soil.Campbell_b + 3.0) / (soil.Campbell_b + 3.0)
    k = soil.Ks * (MFP / soil.CampbellMFP_he)**b3
    return(k)

def dTheta_dH(funcType, soil, H0, H1, z): 
    psi0 = H0 + g*z
    psi1 = H1 + g*z
    if (abs(psi1-psi0) < 1E-5):
        return dTheta_dPsi(funcType, soil, psi0)
    else:
        theta0 = thetaFromPsi(funcType, soil, psi0)
        theta1 = thetaFromPsi(funcType, soil, psi1)
        return (theta1 - theta0) / (psi1 - psi0)
      
def meanK(meanType, k1, k2):
    if (meanType == LOGARITHMIC):
        if (k1 != k2):
            k = (k1-k2) / log(k1/k2)
        else:
            k = k1
    elif (meanType == HARMONIC): 
        k = 2.0 / (1.0 / k1 + 1.0 / k2)
    elif (meanType == GEOMETRIC): 
        k = sqrt(k1 * k2)
    return k

def hydraulicConductivityFromPsi(funcType, soil, psi):
    if (funcType == RESTRICTED_VG):
        psi = abs(psi)
        num = (1. - pow(soil.VG_alpha * psi, soil.VG_m * soil.VG_n)
               *pow(1. + pow(soil.VG_alpha*psi, soil.VG_n), -soil.VG_m))**2
        denom = pow(1. + pow(soil.VG_alpha*psi, soil.VG_n), soil.VG_m * soil.Mualem_L)
        k = soil.Ks * (num / denom)
    elif (funcType == IPPISCH_VG):
        k = NODATA
    elif (funcType == CAMPBELL):
        k = soil.Ks * (soil.Campbell_he / psi)**soil.Campbell_n
    return(k)

########################
###### Vapor-flor ######
########################

def vaporConductivityFromPsiTheta(soil, psi, theta):
    humidity = exp(mw*psi/(R*T)) 
    k = 0.66*(soil.thetaS-theta)*dv*vp*humidity*mw/(R*T)
    return(k)

def dvapor_dPsi(funcType, soil, psi, theta):
    humidity = exp(mw*psi/(R*T)) 
    capacity_vapor = (soil.thetaS-theta)*vp*humidity*(mw/(R*T))-dTheta_dPsi(funcType,soil,psi)*vp*humidity
    return(capacity_vapor)

def vaporFromPsi(soil, psi, theta):
    humidity = exp(mw*psi/(R*T)) 
    vapor = (soil.thetaS - theta)*vp*humidity
    return(vapor)

def evaporation_flux(psi):
    h_s = exp(mw*psi/(R*T)) 
    return(E_p * (h_s-h_a)/(1.-h_a))
	