import math
import numpy as np
from scipy import sparse

# import py_rootbox as rb   
from rb_tools import *
from numpy.linalg.linalg import norm


# Units
#
# M = kg
# L = cm
# T = s
# 
def xylem_flux_ls(seg, nodes, sLen, sType, sRad, kr, kz):
    
    rho = 1e-3 # kg / cm  
    g = 9.8 *100 # cm / s^2 
        
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
    
        n1 = v2v(nodes[i])
        n2 = v2v(nodes[j])
        mid = 0.5*(n1+n2)
        
        v = n2-n1
        v = v / norm(v)
        # v[2] = 0 # TODO
        
        p_s = -60 # TODO psoil(mid) 
        
        a = sRad[c]
        l = sLen[c]
        t = int(sType[c])
        
        bi =  -2.*a*math.pi*l*kr[t]*p_s - kz[t]*rho*g*abs(v[2])# Eqn 9
        cii = -a*math.pi*l*kr[t] - kz[t]/l # Eqn 10
        cij = -a*math.pi*l*kr[t] + kz[t]/l # Eqn 11
        
        # edge ij
        b[i] +=  bi           
        
        I[k] = i
        J[k] = i     
        V[k] += cii
        k += 1                
        
        I[k] = i
        J[k] = j        
        V[k] += cij
        k += 1 
        
        # edge ji
        i,j = j, i
        
        b[i] += bi          

        I[k] = i
        J[k] = i  
        V[k] += cii    
        k += 1                
        
        I[k] = i
        J[k] = j        
        V[k] += cij
        
        k += 1 
        
    Q = sparse.coo_matrix((V,(I,J)))
    Q = sparse.csr_matrix(Q) # Sparse row matrix seems the most reasonable to solve Qx = b iteratively
    
    return (Q, b)


#
#
#
def xylem_flux_bc_dirichlet(Q, b, seg, d):
    c = 0
    for s in seg:
        i = s[0]      
        e0 = np.zeros((1,Q.shape[1])) # build zero vector
        Q[i-1,:] = sparse.csr_matrix(e0) # replace row i with ei
        Q[i-1,i-1] = 1
        b[i-1] = d[c]    
        c += 1

    return Q, b 



