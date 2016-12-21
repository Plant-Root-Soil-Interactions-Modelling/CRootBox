import math
import numpy as np
from scipy import sparse

import py_rootbox as rb   
from rb_tools import *


#
#
#
def xylem_flux(seg, nodes, sLen, sType, sRad, d_):

    kr = [1e-8,1e-8,1e-8,1e-8,1e-8]
    kz = [1e-6,1e-6,1e-6,1e-6,1e-6]
        
    N = len(nodes)
    Ns = len(seg)
    
    I = np.zeros(4*Ns)
    J = np.zeros(4*Ns)    
    V = np.zeros(4*Ns)
    b = np.zeros(N-1)
    
    k = 0

    for c in range(1,Ns): # ignore first segment (bc)
        
        i = seg[c].x-1
        j = seg[c].y-1
        if (i<0):
            print("nooooooo")
        if (j<0):
            print("nooooooo")
        
        
        mid = 0.5*(v2v(nodes[i])+v2v(nodes[j]))
        
        p_soil = -60 # TODO psoil(mid) 
        
        a = sRad[c]
        l = sLen[c]
        t = int(sType[c])
        
        # edge ij
        b[i] -=  2*a*math.pi*l*kr[t]*p_soil # term for gravity is neglected           

        V[k] += -a*math.pi*l*kr[t] - kz[t]/l
        I[k] = i
        J[k] = i        
        k += 1                
        
        V[k] +=  -a*math.pi*l*kr[t] + kz[t]/l;
        I[k] = i
        J[k] = j        
        k += 1 
        
        # edge ji
        i,j = j, i
        b[i] -=  2*a*math.pi*l*kr[t]*p_soil # term for gravity is neglected           

        V[k] += -a*math.pi*l*kr[t] - kz[t]/l
        I[k] = i
        J[k] = i        
        k += 1                
        
        V[k] +=  -a*math.pi*l*kr[t] + kz[t]/l;
        I[k] = i
        J[k] = j        
        k += 1 
        
    Q = sparse.coo_matrix((V,(I,J)))
    
    a = sRad[0] # scaling factor for dir
    l = sLen[0]
    t = int(sType[0])    
    d = -a*math.pi*l*kr[t] + kz[t]/l
    b[0] = b[0] - d*d_;
    
    return (Q, b)

