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
    
    I = np.zeros(Ns);
    J = np.zeros(Ns);    
    V = np.zeros(Ns);
    
    for c in range(0,Ns):
        
        i = seg[c].x
        j = seg[c].y
        I[c] = i
        J[c] = j        
        
        mid = 0.5*(v2v(nodes[i])+v2v(nodes[j]))
        p_soil = -60 # TODO psoil(mid) 
        a = sRad[c]
        l = sLen[c]
        t = sType[c]
        
        # edge ij
        b(i) = b(i) - 2*a*math.pi*l*kr(t)*p_soil # term for gravity is neglected           
        Q(i,i) = Q(i,i) + (-a*math.pi*l*kr(t) - kz(t)/l);
        Q(i,j) = Q(i,j) + (-a*math.pi*l*kr(t) + kz(t)/l);
   
        # edge ji
        i,j = j, i
        b(i) = b(i) - 2*a*pi*l*kr(t)*p_soil; % term for gravity is neglected   
        Q(i,i) = Q(i,i) + (-a*pi*l*kr(t) - kz(t)/l);
        Q(i,j) = Q(i,j) + (-a*pi*l*kr(t) + kz(t)/l);
        

        
#I = array([0,3,1,0])
#>>> J = array([0,3,1,2])
#>>> V = array([4,5,7,9])
#>>> A = sparse.coo_matrix((V,(I,J)),shape=(4,4))        