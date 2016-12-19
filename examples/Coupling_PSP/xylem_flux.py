import py_rootbox as rb    
import numpy as np
from scipy import sparse

#
# Auxiliary functions that should be moved to py_rootbox
#
def v2v(v): # rb.Vector3 to numpy array 
    return np.array([v.x, v.y, v.z])

def vd2a(vd): # rb.std_vector_double_ to numpy array
    N  = len(vd)
    l = np.zeros(N) 
    for i in range(0,N):
        l[i] = vd[i]
    return l

def z2i(z):
    i = int(round((abs(z)/100)*inf.n))  
    return min(max(i,0),inf.n-1) 



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
        
        mid = 0.5*(v2v(nodes[i])+v2v(nodes[j]))
        p_soil = -60 # TODO psoil(mid) 
        a = sRad[c]
        l = sLen[c]
        t = sType[c]
        
        
        b(i) = b(i) - 2*a*math.pi*l*kr(t)*p_soil # term for gravity is neglected           
        Q(i,i) = Q(i,i) + (-a*math.pi*l*kr(t) - kz(t)/l);
        Q(i,j) = Q(i,j) + (-a*math.pi*l*kr(t) + kz(t)/l);
   
   % edge ji
   jj = j;
   j = i;
   i = jj;
   b(i) = b(i) - 2*a*pi*l*kr(t)*p_soil; % term for gravity is neglected   
   Q(i,i) = Q(i,i) + (-a*pi*l*kr(t) - kz(t)/l);
   Q(i,j) = Q(i,j) + (-a*pi*l*kr(t) + kz(t)/l);
        
        I[c] = i
        J[c] = j
        
#I = array([0,3,1,0])
#>>> J = array([0,3,1,2])
#>>> V = array([4,5,7,9])
#>>> A = sparse.coo_matrix((V,(I,J)),shape=(4,4))        