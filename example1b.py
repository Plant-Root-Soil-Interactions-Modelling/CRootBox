import py_rootbox as rb

rootsystem = rb.RootSystem()

# Open plant and root parameter from a file
name = "Anagallis_femina_Leitner_2010" 
rootsystem.openFile(name)

# Create and set geometry

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

# pick 1, 2, or 3
rootsystem.setGeometry(rotatedRhizotron)  # soilcore, rhizotron, or rotatedRhizotron

# Initialize
rootsystem.initialize() # make sure to call rootsystem.setGeometry before initialize()

# Simulate
rootsystem.simulate(60)

# Export final result (as vtp)
rootsystem.write("results/"+name+".vtp")  # roots are exported as polyline

# Export container geometry as Paraview Python script (run file in Paraview by Tools->Python Shell, Run Script)    
rootsystem.write("results/"+name+".py") 
