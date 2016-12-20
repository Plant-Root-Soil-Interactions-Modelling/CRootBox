import py_rootbox as rb   
 
import numpy as np
from scipy import sparse


#
#
#
def xylem_flux(seg, nodes, sLen, sType, sRad):

    kr = [1e-8,1e-8,1e-8,1e-8,1e-8]
    kz = [1e-6,1e-6,1e-6,1e-6,1e-6]
        
    N = len(nodes)
    Ns = len(seg)
    
    I = np.zeros(4*Ns)
    J = np.zeros(4*Ns)    
    V = np.zeros(4*Ns)
    b = np.zeros(N)
    
    k = 0

    for c in range(0,Ns):
        
        i = seg[c].x
        j = seg[c].y 
        
        mid = 0.5*(v2v(nodes[i])+v2v(nodes[j]))
        p_soil = -60 # TODO psoil(mid) 
        a = sRad[c]
        l = sLen[c]
        t = sType[c]
        
        # edge ij
        b[i] -=  2*a*math.pi*l*kr(t)*p_soil # term for gravity is neglected           

        V[k] += -a*math.pi*l*kr[t] - kz[t]/l
        I[k] = i
        J[k] = i        
        k += 1                
        
        V[k] +=  -a*math.pi*l*kr(t) + kz(t)/l;
        I[k] = i
        J[k] = j        
        k += 1 
        
        # edge ji
        i,j = j, i
        b[i] -=  2*a*math.pi*l*kr(t)*p_soil # term for gravity is neglected           

        V[k] += -a*math.pi*l*kr[t] - kz[t]/l
        I[k] = i
        J[k] = i        
        k += 1                
        
        V[k] +=  -a*math.pi*l*kr(t) + kz(t)/l;
        I[k] = i
        J[k] = j        
        k += 1 
        
    Q = sparse.coo_matrix((V,(I,J)))
    return (Q, b)

