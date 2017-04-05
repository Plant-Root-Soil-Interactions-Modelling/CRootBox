#
# Example 5
# 
# Hydrotropism ,
# proof of concept, with a static soil water content
# 
import py_rootbox as rb

rootsystem = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010" 

#
# Plant and root parameter from a file
# 
rootsystem.openFile(name);

# "manually" set tropism to hydrotropism for the first ten root types
for i in range(0,10):    
    rootsystem.getRootTypeParameter(i+1).tropismT = rb.TropismType.hydro;
    rootsystem.getRootTypeParameter(i+1).tropismN = 3; # N
    rootsystem.getRootTypeParameter(i+1).tropismS = 0.4; # sigma

# check if it worked
# for i in range(0,10):
#     print(rootsystem.getRootTypeParameter(i+1))
 
print (rootsystem.getRootSystemParameter().seedPos)
   
#
# Static soil property
#    
sideBox = rb.SDF_PlantBox(10,20,50)
left = rb.SDF_RotateTranslate(sideBox, rb.Vector3d(-4.99,0,0))
right = rb.SDF_RotateTranslate(sideBox, rb.Vector3d(4.99,0,0))
leftright = rb.SDF_Union(left,right)
rootsystem.setGeometry(leftright)  # for vizualisation

maxS = 0.7 # maximal saturation
minS = 0.1 # minimal saturation
slope = 20 # [cm] linear gradient between min and ma
soilprop = rb.SoilPropertySDF(left, maxS, minS, slope)

# set the soil properties before calling initialize
rootsystem.setSoil(soilprop);

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
rootsystem.write("results/"+name + ".vtp", rb.OutputType.polylines);

#
# Export container geometry as Paraview Python script (run file in Paraview by Tools->Python Shell, Run Script)
# 
rootsystem.write("results/"+name + ".py",0);

print("Finished with a total of " + str(rootsystem.getNumberOfNodes()) + " nodes\n")
