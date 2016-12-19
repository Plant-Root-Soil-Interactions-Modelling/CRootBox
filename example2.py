#
# Example 2
#
# Same as Example 1, but with a plant container confining root growth
#
# Example container are:
# 1. Cylindrical soil core,
# 2. Flat square rhizotron,
# 3. Rotated square rhizotron
#
# Additionally, exports the confining geometry as paraview pyhton script
# (run file in Paraview by Tools->Python Shell, Run Script)
#
import py_rootbox as rb

rootsystem = rb.RootSystem()
name = "lupine2014" 

#
# Open plant and root parameter from a file
#
rootsystem.openFile(name)

#
# set geometry
#

# 1.creates a cylindrical soil core with top radius 5 cm, bot radius 5 cm, height 50cm, not square but circular
soilcore = rb.SDF_PlantContainer(5,5,40,False)

# 2. creates a square 27*27 cm containter with height 1.5 cm (used in parametrisation experiment
rhizotron = rb.SDF_PlantBox(1.4,27,27);

# 3. creates a square rhizotron r*r, with height h, rotated around the x-axis for angle alpha
r = 20
h = 4
alpha = 45
rhizotron2 = rb.SDF_PlantContainer(r,r,h,True);
posA = rb.Vector3d(0,r,-h/2); #  origin before rotation
A = rb.Matrix3d.rotX(alpha/180.*3.14)
posA = A.times(posA) # origin after rotation
rotatedRhizotron = rb.SDF_RotateTranslate(rhizotron2,alpha,0,posA.times(-1));

rootsystem.setGeometry(rotatedRhizotron)  # soilcore, rhizotron, or rotatedRhizotron

#
# Initialize
#
rootsystem.initialize()

#
# Simulate
#
simtime = 60  # or 20, 40, 60 days
dt = 1 # try other values here
N = round(simtime/dt) # steps

for i in range(0,int(N)):
    rootsystem.simulate(dt);

#
# Export final result (as vtp)
#
rootsystem.write(name+".vtp",rb.OutputType.polylines) # use ot_polylines for nicer visualization, ot_segments for animations

#
# Export container geometry as Paraview Python script (run file in Paraview by Tools->Python Shell, Run Script)
#    
rootsystem.write(name+".py", 0) 

print("Finished with a total of "+str(rootsystem.getNumberOfNodes())+ " nodes")
