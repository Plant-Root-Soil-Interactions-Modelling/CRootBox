#
# Example 5b
# 
# Cheomtropism ,
# proof of concept, with a static soil water content
#
# soil layer
# 
import py_rootbox as rb

rootsystem = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010" 

#
# Plant and root parameter from a file
# 
rootsystem.openFile(name);

sigma = [0.4, 1., 1., 1., 1. ] * 2
print(sigma)

# "manually" set tropism to hydrotropism for the first ten root types
for i in range(0,10):  
    rootsystem.getRootTypeParameter(i+1).dx = 0.25; # they seem to coarse  
    rootsystem.getRootTypeParameter(i+1).tropismT = rb.TropismType.hydro;
    rootsystem.getRootTypeParameter(i+1).tropismN = 3; # N
    rootsystem.getRootTypeParameter(i+1).tropismS = sigma[i]; # sigma

# check if it worked
# for i in range(0,10):
#     print(rootsystem.getRootTypeParameter(i+1))
 
print (rootsystem.getRootSystemParameter().seedPos)
   
#
# Static soil property
#    
box = rb.SDF_PlantBox(100,100,2)
layer = rb.SDF_RotateTranslate(box, rb.Vector3d(0,0,-16))

maxS = 0.7 # maximal 
minS = 0.1 # minimal 
slope = 5 # [cm] linear gradient between min and ma
soilprop = rb.SoilPropertySDF(layer, maxS, minS, slope)

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
rootsystem.setGeometry(layer)  # for vizualisation
rootsystem.write("results/"+name + ".py",0);

print("Finished with a total of " + str(rootsystem.getNumberOfNodes()) + " nodes\n")
