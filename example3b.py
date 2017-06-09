import py_rootbox as rb
import numpy as np
import matplotlib.pyplot as plt

rootsystem = rb.RootSystem()
name = "Brassica_oleracea_Vansteenkiste_2014" 
rootsystem.openFile(name)
rootsystem.initialize() 
rootsystem.simulate(120, True) 
rootsystem.write(name+".vtp")

#
# Soil core analysis
#
r = 10
depth = 100.
layers = 100
soilcolumn = rb.SDF_PlantContainer(r,r,depth,False) # in the center of the root
soilcolumn2 = rb.SDF_RotateTranslate(soilcolumn,0,0,rb.Vector3d(10,0,0)) # shift 10 cm
geom = soilcolumn2

z_=np.linspace(0,-1*depth,layers)
fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(16,8))
for a in axes:    
    a.set_xlabel('RLD (cm/cm)')
    a.set_ylabel('Depth (cm)')

# make a root length distribution
ana = rb.SegmentAnalyser(rootsystem)
rl_ = ana.distribution(rb.ScalarType.length,0.,depth,layers,True)
axes[0].set_title('All roots (120 days)')
axes[0].plot(rl_,z_)

# make a root length distribution along the soil core
ana = rb.SegmentAnalyser(rootsystem)
ana.crop(geom)
ana.pack()
rl_ = ana.distribution(rb.ScalarType.length,0.,depth,layers,True)
axes[1].set_title('Soil core (120 days)')
axes[1].plot(rl_,z_)

# how it looked after 30 days?
ana = rb.SegmentAnalyser(rootsystem)
ana.filter(rb.ScalarType.time,0,30)
ana.crop(geom)
ana.pack()
rl_ = ana.distribution(rb.ScalarType.length,0.,depth,layers,True)
axes[2].set_title('Soil core (30 days)')
axes[2].plot(rl_,z_)

# only laterals?
ana = rb.SegmentAnalyser(rootsystem)
ana.filter(rb.ScalarType.type,2) # assuming laterals are of type 2
ana.crop(geom)
ana.pack()
rl_ = ana.distribution(rb.ScalarType.length,0.,depth,layers,True)
axes[3].set_title('Soil core, lateral roots (120 days)')
axes[3].plot(rl_,z_)

fig.subplots_adjust()
plt.savefig("results/example_3b.png")
plt.show()
