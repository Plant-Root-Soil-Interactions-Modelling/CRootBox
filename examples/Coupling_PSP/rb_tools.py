import numpy as np



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

def vv2a(vd): # rb.std_vector_Vector3_ to numpy array
    N  = len(vd)
    l = np.zeros((N,3)) 
    for i in range(0,N):
        l[i,:] = [vd[i].x,vd[i].y,vd[i].z]
    return l


def z2i(z,n): # maps z to equidistant mesh
    i = int(round((abs(z)/100)*n))  
    return min(max(i,0),n-1) 

def plotRSinit(ax):
    ax.clear()
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_zlim(-1, 0)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")     

def plotRS(ax,seg,nodes): # plots the root system as line plot (too slow)
    plotRSinit(ax)
    scale = 0.01;
    for s in seg:
        n1 = v2v(nodes[s.x])*scale
        n2 = v2v(nodes[s.y])*scale
        ax.plot((n1[0],n2[0]),(n1[1],n2[1]),(n1[2],n2[2]),'k-')
        
def plotRSscatter(ax,nodes): # plots the root system nodes (rather slow)
    plotRSinit(ax)
    scale = 0.01
    n = vv2a(nodes)
    ax.scatter(n[:,0]*scale,n[:,1]*scale,n[:,2]*scale)