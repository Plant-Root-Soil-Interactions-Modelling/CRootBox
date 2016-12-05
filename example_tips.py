import py_rootbox as rb

import numpy as np
import matplotlib.pyplot as plt

from math import sqrt



rootsystem = rb.RootSystem()
name = "zeamays_test" # anagallis2010, zeamays_test

#
# Open plant and root parameter from a file
#
rootsystem.openFile(name,"modelparameter/")
# rootsystem.writeParameters() # not exposed to python yet

#
# Initialize
#
rootsystem.initialize(4,5) # TODO expose default values

#
# Simulate
#
simtime = 120;
rootsystem.simulate(simtime);
    
#
# Analyse root tips
# 
tips = rootsystem.getRootTips()
notips = sum(1 for _ in tips) # is there a more clever way?

print('Number of root tips is '+str(notips))

tipcoords = np.zeros((notips,3)) 
z_ = np.zeros(notips)
l_ = np.zeros(notips)

#
# Copy everything we need
#
c=0;
for t in tips:
    tipcoords[c,0] = t.x
    tipcoords[c,1] = t.y
    tipcoords[c,2] = t.z  
    z_[c] = t.z
    l_[c] = sqrt(t.x*t.x + t.y*t.y)  
    c+=1
#
# Figure params
#
lbins = 50
lrange = (0,50)
zbins = 120
zrange = (-120,0)

#
# Plot histograms
#
plt.figure(1)
plt.hist(l_, bins=lbins, range=lrange)  
plt.title("Root tip radial distance")
plt.show(False)
hl, bins = np.histogram(l_, bins=lbins, range=lrange)
print(hl);

np.savetxt('radialdistribution.txt',hl); # save results

plt.figure(2)
plt.hist(z_, bins=zbins, range=zrange)  
plt.title("Root tip depth")
plt.show()
hz, bins = np.histogram(z_, bins=zbins, range=zrange)
print (hz);

np.savetxt('depthdistribution.txt',hz); # save results

