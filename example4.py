#
# Example 4
# 
# More complex geometries
# 
# 1. a split pot experiment
# 2. rhizotubes as obstacles (from Example_Rhizotubes.m)
# 
# Additionally, exports the confining geometry as paraview pyhton script
# (run file in Paraview by Tools->Python Shell, Run Script)
# 
import py_rootbox as rb

rootsystem = rb.RootSystem()
name = "maize" 

#
# Plant and root parameter from a file
#
rootsystem.openFile(name,"modelparameter/");

#
# 1. A split pot experiment
#
topBox = rb.SDF_PlantBox(22,20,5)
sideBox = rb.SDF_PlantBox(10,20,35)

left = rb.SDF_RotateTranslate(sideBox, rb.Vector3d(-6,0,-5))
right = rb.SDF_RotateTranslate(sideBox, rb.Vector3d(6,0,-5))

box_ = rb.std_vector_SDF_()
box_.append(topBox)
box_.append(left)
box_.append(right)
splitBox = rb.SDF_Union(box_);

#
# 2. Rhizotubes as obstacles 
#

# Box
boxX = 96
boxY = 126
boxZ = 130
box = rb.SDF_PlantBox(boxX,boxY,boxZ);

# A single Rhizotube
r = 2*3.2; # cm
rhizotube = rb.SDF_PlantContainer(r,r,96,False);
rhizoX = rb.SDF_RotateTranslate(rhizotube, 90, rb.SDF_Axis.yaxis, rb.Vector3d(boxX/2,0,0));

# The experimental setting
rhizotubes_ = rb.std_vector_SDF_()
y_ = ( -30, -18, -6, 6, 18, 30 )
z_ = ( -10, -20, -40, -60, -80, -120 )

for i in range(0,len(y_)):
    tube = rb.SDF_RotateTranslate(rhizoX, rb.Vector3d(0,y_[i],z_[i]))
    rhizotubes_.append(tube)
    
# Final geometry    
rhizotubes = rb.SDF_Union(rhizotubes_)
rhizoTube = rb.SDF_Difference(box, rhizotubes)

rootsystem.setGeometry(splitBox); # splitBox, or rhizoTube

#
# Initialize
#
rootsystem.initialize(4,5); # it is important to call initialize() after setGeometry()

#
# Simulate
# 
simtime = 100; # e.g. 30 or 60 days
dt = 1;
N = round(simtime/dt);
for _ in range(0,N):        
    rootsystem.simulate(dt)   

#     
# Export results (as vtp)
#    
rootsystem.write(name + ".vtp", rb.OutputType.polylines);

#
# Export container geometry as Paraview Python script (run file in Paraview by Tools->Python Shell, Run Script)
# 
rootsystem.write(name + ".py",0);

print("Finished with a total of " + str(rootsystem.getNumberOfNodes()) + " nodes\n")
