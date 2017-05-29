#
# Compares numerical approximations with altering resolutions dt, dx (L86 - )
#
# plots root tip radial distance
# plots root tip depth 
# x-y tip density map  
# 
# computes in parallel to enable a lot of runs
#
import py_rootbox as rb

from multiprocessing import Pool
from itertools import product

from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats



def v2a(vd): # rb.std_vector_double_ to numpy array    
    l = np.zeros(len(vd)) 
    for i in range(0,len(vd)):
        l[i] = vd[i]
    return l

def a2v(a): #  numpy array to rb.std_vector_double
    l = rb.std_vector_double_()
    for d in a:
        l.append(d)
    return l

def a2i(a): #  numpy array to rb.std_int_double
    l = rb.std_vector_int_()
    for i in a:
        l.append(i)
    return l

def vv2a(vd): # rb.std_vector_Vector3_ to numpy array
    N  = len(vd)
    l = np.zeros((N,3)) 
    for i in range(0,N):
        l[i,:] = [vd[i].x,vd[i].y,vd[i].z]
    return l

def simOnce(name, simtime, lbins, lrange, zbins, zrange, dx, dt):
    # simulation
    rs = rb.RootSystem()
    rs.openFile(name,"modelparameter/")
    for i in range(0,10):    
        rs.getRootTypeParameter(i+1).dx = dx    
    rs.initialize() 
    N = round(simtime/dt)
    print("simOnce")
    for i in range(0,N):
        rs.simulate(dt,True);
    # analysis
    img = np.zeros((2*lbins,2*lbins))
    nodes = vv2a(rs.getNodes());
    tips = rs.getRootTips()
    notips = len(tips)
    z_ = np.zeros(notips)
    l_ = np.zeros(notips)
    c=0;
    for t in tips:
        x = nodes[t,0]
        y = nodes[t,1]
        z = nodes[t,2]  
        # tip top view
        i = np.around((x/lrange[1])*lbins+lbins)
        j = np.around((y/lrange[1])*lbins+lbins)
        i = min(i,2*lbins-1)
        j = min(j,2*lbins-1)        
        i = int(max(i,0))
        j = int(max(j,0))        
        img[i,j] += 1.
        # depth, and length distribution
        z_[c] = z
        l_[c] = sqrt(x*x + y*y)          
        c+=1   
        hl, bins = np.histogram(l_, bins=lbins, range=lrange)
        hz, bins = np.histogram(z_, bins=zbins, range=zrange)        
    return hl,hz,img

# Params
dx = 2
dt = 60

runs = 300
# name = "Lupinus_albus_Leitner_2014"
# name = "Zea_mays_1_Leitner_2010"
name = "wheat"
simtime = 240

# Histogram params
lbins = 40
lrange = (0,20)
zbins = 120
zrange = (-120,0)


pool = Pool() #defaults to number of available CPU's
chunksize = 20 #this may take some guessing ... take a look at the docs to decide
output =  pool.starmap(simOnce,[(name,simtime,lbins,lrange,zbins,zrange,dx,dt)]*runs)

allL, allZ,tiptop = output[0];
for i in range(1,runs):
    L, Z, img  = output[i];
    allL = np.vstack((allL,L))
    allZ = np.vstack((allZ,Z))
    tiptop += img

meanZ = np.mean(allZ,0)
semZ = stats.sem(allZ,0)
meanL = np.mean(allL,0)
semL = stats.sem(allL,0)

plt.figure(1)
plt.errorbar(np.linspace(lrange[0],lrange[1],lbins), meanL, semL, linestyle='None', marker='^')
plt.title("Root tip radial distance")
plt.show(False)

plt.figure(2)
plt.errorbar(np.linspace(zrange[0],zrange[1],zbins), meanZ, semZ, linestyle='None', marker='^')
plt.title("Root tip depth")

plt.figure(3)
plt.title("Tip top view")
plt.imshow(tiptop)

plt.show()



