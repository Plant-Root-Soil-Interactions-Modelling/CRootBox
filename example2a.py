import py_rootbox as rb
import math

rootsystem = rb.RootSystem()
name = "Zea_mays_4_Leitner_2014"

# Open plant and root parameter from a file
rootsystem.openFile(name)

# 1. creates a square rhizotron r*r, with height h, rotated around the x-axis for angle alpha
r = 20
h = 4
alpha = 45
rhizotron2 = rb.SDF_PlantContainer(r,r,h,True)
posA = rb.Vector3d(0,r,-h/2) # origin before rotation
A = rb.Matrix3d.rotX(alpha/180.*math.pi)
posA = A.times(posA) # origin after rotation
rotatedRhizotron = rb.SDF_RotateTranslate(rhizotron2,alpha,0,posA.times(-1))

# 2. A split pot experiment
topBox = rb.SDF_PlantBox(22,20,5)
sideBox = rb.SDF_PlantBox(10,20,35)
left = rb.SDF_RotateTranslate(sideBox, rb.Vector3d(-6,0,-5))
right = rb.SDF_RotateTranslate(sideBox, rb.Vector3d(6,0,-5))
box_ = rb.std_vector_SDF_()
box_.append(topBox)
box_.append(left)
box_.append(right)
splitBox = rb.SDF_Union(box_)

# 3. Rhizotubes as obstacles 
# Box
box = rb.SDF_PlantBox(96,126,130)
# A single Rhizotube
r = 2*3.2 # cm
rhizotube = rb.SDF_PlantContainer(r,r,96,False)
rhizoX = rb.SDF_RotateTranslate(rhizotube, 90, rb.SDF_Axis.yaxis, rb.Vector3d(96/2,0,0))
# The experimental setting
rhizotubes_ = rb.std_vector_SDF_()
y_ = ( -30, -18, -6, 6, 18, 30 )
z_ = ( -10, -20, -40, -60, -80, -120 )
tube = [] 
for i in range(0,len(y_)):
    v = rb.Vector3d(0,y_[i],z_[i])
    tube.append(rb.SDF_RotateTranslate(rhizoX, v))
    rhizotubes_.append(tube[i])    
# Final geometry    
rhizotubes = rb.SDF_Union(rhizotubes_)
rhizoTube = rb.SDF_Difference(box, rhizotubes)

# set geometry
rootsystem.setGeometry(rotatedRhizotron) # rotatedRhizotron, splitBox, or rhizoTube

# Initialize
rootsystem.initialize() 

# Simulate
rootsystem.simulate(60) # days   

# Export results (as vtp)    
rootsystem.write("results/example_2a.vtp")

# Export container geometry as Paraview Python script 
rootsystem.write("results/example_2a.py")

