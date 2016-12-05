import py_rootbox as rb

import numpy as np
import matplotlib.pyplot as plt

from math import sqrt


output = np.zeros((N,N))
pool = Pool() #defaults to number of available CPU's
chunksize = 20 #this may take some guessing ... take a look at the docs to decide
for ind, res in enumerate(pool.imap(Fun, product(xrange(N), xrange(N))), chunksize):
    output.flat[ind] = res

def simOnce(name,simtime,lbins,lrange,zbins,zrange):
    rootsystem = rb.RootSystem()
    rootsystem.openFile(name,"modelparameter/")
    rootsystem.initialize(4,5) # TODO expose default values
    rootsystem.simulate(simtime);
    tips = rootsystem.getRootTips()
    notips = sum(1 for _ in tips) # is there a more clever way?
    tipcoords = np.zeros((notips,3)) 
    z_ = np.zeros(notips)
    l_ = np.zeros(notips)
    c=0;
    for t in tips:
        tipcoords[c,0] = t.x
        tipcoords[c,1] = t.y
        tipcoords[c,2] = t.z  
        z_[c] = t.z
        l_[c] = sqrt(t.x*t.x + t.y*t.y)          
        c+=1   
        hl, bins = np.histogram(l_, bins=lbins, range=lrange)
        hz, bins = np.histogram(z_, bins=zbins, range=zrange)
    return hl,hz
    

# Params
runs = 100
name = 'zeamays_test'
simtime = 120

# Histogram params
lbins = 50
lrange = (0,50)
zbins = 120
zrange = (-120,0)

allL, allZ = simOnce(name,simtime,lbins,lrange,zbins,zrange);
# print(np.shape(allZ)) # size is shape
for i in range(1,runs):
    L, Z  = simOnce(name,simtime,lbins,lrange,zbins,zrange);
    allL = np.vstack((allL,L));
    allZ = np.vstack((allZ,Z));
    print(np.shape(allL)) # size is shape    
    print(np.shape(allZ)) # size is shape
        
meanZ = np.mean(allZ,0)
stdZ = np.std(allZ,0)
meanL = np.mean(allL,0)
stdL = np.mean(allL,0)

print(np.shape(meanL))
print(np.shape(stdL))

plt.figure(1)
plt.errorbar(np.linspace(lrange[0],lrange[1],lbins), meanL, stdL, linestyle='None', marker='^')
plt.title("Root tip radial distance")
plt.show(False)

plt.figure(2)
plt.errorbar(np.linspace(zrange[0],zrange[1],zbins), meanZ, stdZ, linestyle='None', marker='^')
plt.title("Root tip depth")
plt.show()

np.savetxt('radialdistribution.txt',meanL); # save results
np.savetxt('depthdistribution.txt',meanZ); # save results



